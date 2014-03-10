% test for farthest point sampling

n = 300;

name = 'bump';
name = 'map';
name = 'stephanodiscusniagarae';
name = 'cavern';
name = 'gaussian';
name = 'road2';
name = 'binary';
name = 'constant';
name = 'mountain';
name = 'anisotropic';

mask_name = 'disk';
mask_name = '';
mask = ones(n);
if not(isempty(mask_name))
    options.radius = .48;
    mask = rescale( load_image(mask_name,n,options) );
end
options.constraint_map = mask;
options.constraint_map(mask>.5) = +Inf;
options.constraint_map(mask<=.5) = -Inf;


rep = ['results/farthest-sampling/' name '/'];
if exist(rep)~=7
    mkdir(rep);
end

if strcmp(name, 'anisotropic')
    % anisotropic tensor field
    aniso_type = 'fixed';
    aniso_type = 'varying';
    if strcmp(aniso_type, 'fixed')
        u = [1 .5];
        U = repmat( reshape(u,1,1,2),[n n 1] );
    elseif strcmp(aniso_type, 'varying')
        U = randn(n,n,2);
        sigma = 40*(n/200);
        for it=1:10
            U = perform_vf_normalization( perform_blurring(U, sigma) );
        end
    end
    U = perform_vf_normalization( U );
    % 3D field
    aniso = .01;
    V = cat(3, -U(:,:,2), U(:,:,1)); % orthogonal vector
    W = perform_tensor_recomp(U,V, ones(n),ones(n)*aniso );    
    if exist('perform_lic')
        % compute a cool 2D lic texture
        options.isoriented = 0;
        M0 = perform_blurring(randn(n),0);
        M0 = perform_histogram_equalization( M0, 'linear');
        options.histogram = 'linear';
        options.dt = 0.4; options.M0 = M0;
        options.verb = 1; options.flow_correction = 1;
        options.niter_lic = 2;
        w = 15;
        w = 30;
        M = perform_lic(U, w, options);
    else
        M = U; % background image
    end
    clf;
    plot_tensor_field(W,M);
    saveas(gcf, [rep name '-tensors.png'], 'png');
elseif strcmp(name, 'binary')
    W = rescale(W,0.5,1); M = W;
    M(1) = min(M(:))-.3; M(2) = max(M(:))+.3;
else
    [M,W] = load_potential_map(name, n);
end

warning off;
imwrite(rescale(M), [rep name '-metric.png'], 'png');
warning on;


% plot sampling location
ms = 12; lw = 1.5;
i = 0;

I = find(mask==1);
[vx vy] = ind2sub(size(W), I(1));
landmark = [vx;vy];


for nbr_landmarks = [100 200 500 800 1200] % 10 20 40 
    i = i+1;
    disp('Perform farthest point sampling');
    landmark = perform_farthest_point_sampling( W, landmark, nbr_landmarks-size(landmark,2), options );
    
    % compute the associated triangulation
    [D,Z,Q] = perform_fast_marching(W, landmark);

    Q1 = zeros(n+2)-1; Q1(2:end-1,2:end-1) = Q;                    
    faces = compute_voronoi_triangulation(Q1,landmark);
    edges = compute_edges(faces);
    I = find(edges(1,:)>0 & edges(2,:)>0 );
    edges = edges(:,I);
    
    sel = randperm(nbr_landmarks);  % shuffle colors
    Q = sel(Q);
    
    clf; 
    hold on;
    imageplot(Q'); 
    plot(landmark(1,:), landmark(2,:), 'r.', 'MarkerSize', ms);
    colormap jet(256);
    hold off;
    saveas(gcf, [rep name '-voronoi-' num2string_fixeddigit(nbr_landmarks,3) '.png'], 'png');
    
    % display sampling + distance
    D = perform_histogram_equalization(D,linspace(0,1,n^2));   
    clf;
    hold on;
    imageplot(D');
    plot(landmark(1,:), landmark(2,:), 'r.', 'MarkerSize', ms);
    hold off;
    colormap jet(256);
    saveas(gcf, [rep name '-sampling-' num2string_fixeddigit(nbr_landmarks,3) '.png'], 'png');
    
    % display triangulation
    clf;
    hold on;
    if strcmp(name, 'anisotropic') && ~exist('perform_lic')
        plot_tensor_field(W, M);
    else
        imageplot( M');
    end
    if not(isempty(edges))
        h = plot_edges(edges, landmark, 'b');
        set(h, 'LineWidth',lw);
    end
    plot(landmark(1,:), landmark(2,:), 'r.', 'MarkerSize', ms);
    hold off;
    axis tight; axis image; axis off;
    colormap gray(256);
    saveas(gcf, [rep name '-triangulation-' num2string_fixeddigit(nbr_landmarks,3) '.png'], 'png');
end