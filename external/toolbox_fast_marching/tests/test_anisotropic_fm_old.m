% test for a simple anisotropic metric
name = 'fixed-2d';
name = 'fixed-3d';
name = 'varying-2d';

rep = 'results/anisotropic-fm/';
if not(exist(rep))
    mkdir(rep);
end

% create the main direction of the field
switch name
    case 'fixed-3d' % Fixed 3D tensor field        
        n = 20;
        s = [n n n]; % size
        % main direction of the tensor
        u = [1 .1 .1];
        U = repmat( reshape(u,1,1,1,3),[n n n 1] );
    case 'fixed-2d' % Fixed 2D tensor field        
        n = 40;
        s = [n n 1];
        % main direction of the tensor
        u = [1 .1];
        U = repmat( reshape(u,1,1,2),[n n 1] );
    case 'varying-2d' % spacially varying 2D field
        n = 200;
        % create 2D vector field
        s = [n n 1];        
        U = randn(n,n,2);
        sigma = 30;
        for it=1:10
            U = perform_vf_normalization( perform_blurring(U, sigma) );
        end
        
end
U = perform_vf_normalization( U );


% test for various degree of anisotropy
aniso_list = [0.01 0.05 .1 .2 .5 1];

nstart = 1;
nstart = 8;
if nstart==1
    start_points = [n/2 n/2];
else
    start_points = rand(2,nstart);
    start_points(1,:) = rescale(start_points(1,:),.1,.9);
    start_points(2,:) = rescale(start_points(2,:),.1,.9);
    start_points = round( start_points*(n-1)+1 );
end
if strcmp(name(end-1:end), '3d')
    start_points(end/2) = ceil(s(end)/2); 
end

for ianiso = 1:length(aniso_list)
    
    aniso = aniso_list(ianiso);

    % use cross product to compute the 2 remaining orthogonal directions
    if strcmp(name(end-1:end), '2d')
        % 3D field
        V = cat(3, -U(:,:,2), U(:,:,1)); % orthogonal vector
        T = perform_tensor_recomp(U,V, ones(n),ones(n)*aniso );
    else
        % 3D field
        U = cat(5, U, randn(s(1),s(2),s(3),3,2));
        U(:,:,:,:,3) = cross( U(:,:,:,:,1),U(:,:,:,:,2), 4 );
        U(:,:,:,:,3) = perform_vf_normalization( U(:,:,:,:,3) );
        U(:,:,:,:,2) = cross( U(:,:,:,:,1),U(:,:,:,:,3), 4 );
        U(:,:,:,:,2) = perform_vf_normalization( U(:,:,:,:,2) );
        Lambda = ones(s(1),s(2),s(3),3);
        Lambda(:,:,:,2:3) = aniso;
        T = perform_tensor_decomp_3d(U,Lambda);
    end

    D = perform_fast_marching(T, start_points);

    D1 = D(:,:,ceil(s(end)/2));
    if strcmp(name(end-1:end),'3d')
        T1 = T(:,:,ceil(s(end)/2), 1:2,1:2);
    else
        T1 = T;
    end
        
    clf;
    imageplot(D1);
    colormap jet(256);
    
    clf;
    hold on;
    options.sub = round(n/15);
    options.color = 'k';
    plot_tensor_field(T1, D1, options);
    h = plot(start_points(2,:),start_points(1,:), 'r.');    
    set(h, 'MarkerSize', 20);
    colormap jet(256);
    
    saveas(gcf, [rep name '-aniso-' num2str(ianiso) '.png'], 'png');

end