% test for a simple anisotropic metric
%name = 'fixed-2d';
%name = 'fixed-3d';
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
        randn('seed', 12345);
        U = randn(n,n,2);
        sigma = (n/200)*30;
        options.bound = 'per';
        for it=1:10
            U = perform_vf_normalization( perform_blurring(U, sigma,options) );
        end
        
end
U = perform_vf_normalization( U );

if exist('perform_lic') && size(U,4)==1
    % compute a cool 2D lic texture
    options.isoriented = 0;
    M0 = perform_blurring(randn(n),0);
    M0 = perform_histogram_equalization( M0, 'linear');
    options.histogram = 'linear';
    options.dt = 0.8; options.M0 = M0;
    options.verb = 1; options.flow_correction = 1;
    options.niter_lic = 1;
    w = 15;
    w = 30;
    M = perform_lic(U, w, options);
    warning off;
    imwrite( rescale(M), [rep name '-texture.png'], 'png' );
    warning on;
else
    M = U; % background image
end

% test for various degree of anisotropy
aniso_list = [.01 .05 .1 .2 .5 1];


%% test for progressive propagation
aniso = .01;
V = cat(3, -U(:,:,2), U(:,:,1)); % orthogonal vector
T = perform_tensor_recomp(U,V, ones(n),ones(n)*aniso );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% one starting point and many geodesics on the boundary
start_points = [n;n]/2;
% end points on the boundary
t = 1:n;
x = [t t*0+n t(end-1:-1:1) t*0+1];
y = [t*0+1 t t*0+n t(end-1:-1:1)];
npaths = 50;
s = round(linspace( 1,length(x), npaths+1) ); s(end) = [];
end_points = cat(1, x(s),y(s));

ms = 30; lw = 3; % display params
for ianiso = 1:length(aniso_list)
    
    % build the tensor field
    aniso = aniso_list(ianiso);
    V = cat(3, -U(:,:,2), U(:,:,1)); % orthogonal vector
    T = perform_tensor_recomp(U,V, ones(n),ones(n)*aniso );
    % propagation
    [D,S,Q] = perform_fast_marching(T, start_points);
    % for sexy display
    D1 = perform_histogram_equalization(D, linspace(0,1,n^2));
    % extract tons of geodesics
    for i=1:npaths
        paths{i} = compute_geodesic(D,end_points(:,i), options);
    end
    % display
    clf; hold on;
    imageplot(D1); axis image; axis off; colormap jet(256);
    for i=1:npaths
        end_point = end_points(:,i);
        h = plot( paths{i}(2,:), paths{i}(1,:), 'k' );
        set(h, 'LineWidth', lw);
        h = plot(end_point(2),end_point(1), '.b');
        set(h, 'MarkerSize', ms);
    end
    h = plot(start_points(2),start_points(1), '.r');
    set(h, 'MarkerSize', ms);
    hold off;
    colormap jet(256);
    axis ij;
    saveas(gcf, [rep name '-geodesics-' num2str(ianiso) '.png'], 'png');

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% many starting point and Voronoi diagrams
p = 3;
x = n/(p*2):n/p:n-n/(p*2);
[Y,X] = meshgrid(x,x);
start_points = round( start_points*(n-1)+1 );
start_points = cat(1, X(:)', Y(:)');    

%% progressive propagation
[D,S,Q] = perform_fast_marching(T, start_points );
dmax_list = linspace( 0,max(D(:)), 15 ); 
dmax_list(1) = []; % dmax_list(end) = Inf;
for i=1:length(dmax_list)
    options.dmax = dmax_list(i);
    [D,S,Q] = perform_fast_marching(T, start_points, options);
    
    I = find(D~=Inf); J = find(D==Inf);
    D(I) = perform_histogram_equalization(D(I), linspace(0,1,length(I)));
    D(J) = 0;
    A = apply_colormap(D, 'jet');
    A(J) = M(J); A(J+n^2) = M(J); A(J+2*n^2) = M(J);
    
    warning off;
    imwrite( A, [rep name '-propagation-' num2str(i)  '.png'], 'png' );
    warning on;
    imageplot(A); drawnow;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test for various anisotropy
nstart = 1;
nstart = 40;
if nstart==1
    start_points = [n/2;n/2];
else
    start_points = rand(2,nstart);
    start_points(1,:) = rescale(start_points(1,:),.1,.9);
    start_points(2,:) = rescale(start_points(2,:),.1,.9);
    p = 8;
    x = n/(p*2):n/p:n-n/(p*2);
    [Y,X] = meshgrid(x,x);
    start_points = round( start_points*(n-1)+1 );
    start_points = cat(1, X(:)', Y(:)');
    nstart = size(start_points,2);
    start_points = start_points(:,randperm(nstart));
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

    [D,S,Q] = perform_fast_marching(T, start_points);

    D1 = D(:,:,ceil(s(end)/2));
    if strcmp(name(end-1:end),'3d')
        T1 = T(:,:,ceil(s(end)/2), 1:2,1:2);
    else
        T1 = T;
    end
        
    D1 = perform_histogram_equalization(D1, linspace(0,1,n^2));
    
    warning off;
    imwrite( apply_colormap(Q, 'jet'), [rep name '-aniso-' num2str(ianiso) '-voronoi.png'], 'png' );
    imwrite( apply_colormap(D1, 'jet'), [rep name '-aniso-' num2str(ianiso) '-distance.png'], 'png' );
    warning on;
    
    clf;
    plot_tensor_field(T1,M);
    saveas(gcf, [rep name '-aniso-' num2str(ianiso) '-tensors.png'], 'png');
    
    clf;
    %hold on;
    options.sub = round(n/15);
%    options.color = 'k';
%    plot_tensor_field(T1, D1, options);
    imageplot(D1);
    h = plot(start_points(2,:),start_points(1,:), 'w.');    
    set(h, 'MarkerSize', 25);
    colormap jet(256);
    
   % saveas(gcf, [rep name '-aniso-' num2str(ianiso) '.png'], 'png');

end