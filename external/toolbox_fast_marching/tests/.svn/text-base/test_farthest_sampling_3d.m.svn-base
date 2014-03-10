% test for farthest point sampling for volumetric mesh

n = 150;

name = 'constant';
name = 'central';
name = 'sphere';
name = 'medical';

rep = ['results/farthest-sampling-3d/' name '/'];
if exist(rep)~=7
    mkdir(rep);
end

switch name
    case 'constant'
        W = ones(n,n,n);
    case {'central' 'sphere'}
        x = linspace(-1,1,n);
        [X,Y,Z] = meshgrid(x, x, x);
        W = sqrt( X.^2+Y.^2+Z.^2 );
        if strcmp(name, 'sphere')
            r = .6;
            W = rescale( -(W<=r) );
        end
        W = rescale(W,.1,1);
    case 'medical'        
        load brain1-crop-256.mat
        M = crop(M,n);
        W = rescale(-abs(M-median(M(:))));
        W = rescale( clamp(W,.6,1),.2,1);
    otherwise
        error('Unknown potential');
end

method = 'geodesic';
method = 'delaunay';

% padd to treat correctly boundary
if strcmp(method, 'geodesic')
    k = 10;
    W1 = perform_image_extension(W,n+2*k, '2side');
end

% plot sampling location
ms = 20; lw = 3;
i = 0;
for npoints = [400 800 1600 3000 5000 10000]
    i = i+1;
    disp('Perform farthest point sampling');
    if i==1
        [x,y,z] = meshgrid([1 n], [1 n], [1 n]);
        vertex = cat(1, x(:)', y(:)', z(:)');
    end
    vertex = perform_farthest_point_sampling( W, vertex, npoints-size(vertex,2) );
    
    % display sampling
    clf;
    hold on;
    h = imageplot(W); 
    view(60,30); zoom(.8);  camlight;
    alphamap(linspace(0,.1,100)); vol3d(h);
    h = plot3(vertex(1,:),vertex(2,:),vertex(3,:), '.');
    set(h, 'MarkerSize', ms);
    hold off;
    saveas(gcf, [rep name '-sampling-' num2string_fixeddigit(npoints,5) '.png'], 'png');
    
    % compute the associated triangulation
    if strcmp(method, 'delaunay')
        faces = delaunay3(vertex(1,:),vertex(2,:),vertex(3,:))';
    else
        vertex1 = vertex+k;
        [D,Z,Q] = perform_fast_marching(W1, vertex1);
        faces = compute_voronoi_triangulation(Q,vertex);
    end
    

    w = [1 1 1]; w = w(:)/sqrt(sum(w.^2));
    t = sum(vertex.*repmat(w,[1 size(vertex,2)]));
    delta = max(t)-min(t);
    offlist = linspace( min(t(:)) + .3*delta, max(t(:)) - .3*delta, 5 );
    for i=1:length(offlist)        
        options.cutting_plane = w;
        options.cutting_offs = offlist(i);
        clf;
        hold on;
        %h = imageplot(W);
        %alphamap(linspace(0,.1,100)); 
        plot_mesh(vertex,faces, options);
        % set(h{2}, 'Marker', 'r.');
        view(60,30); zoom(.8);  camlight;        
        % vol3d(h);
        hold off;
        saveas(gcf, [rep name '-sampling-' num2string_fixeddigit(npoints,5) '-mesh-' num2str(i) '.png'], 'png');
    end
    
end