% test for precise Voronoi cell extraction on meshes, together with Lloyd relaxation

clear options;
path(path, 'toolbox/');
path(path, '../toolbox_graph_data/');
path(path, '../toolbox_graph_data/off/');

% test lloyd algorithm or progressive seeding
test = 'voronoi';
test = 'lloyd';

rep = ['results/' test '-mesh/'];
if not(exist(rep))
    mkdir(rep);
end


name = 'david50kf';
name = 'hand';
name = 'david_head';
name = 'bunny';
name = 'elephant-50kv';

options.name = name;
[vertex,faces] = read_mesh(name);

n = size(vertex,2);


if strcmp(test, 'lloyd')
    nstart = 100;
    niter = 20;
else
    nstart = 5;
    niter = 8;    
    nadd = 15;
    nmax = nstart + (niter-1)*nadd;
end
start_points = round(rand(nstart,1)*(n-1)+1);

cm = rand(nmax,3);
cm = rgb2hsv(cm);
cm(:,2) = 1-(1-cm(:,2))/2;
cm = hsv2rgb(cm);

for i=1:niter
    [Q,DQ, ve, edges_id, lambda] = compute_voronoi_mesh(vertex,faces, start_points, options);
    % display the sampling
    options.voronoi_edges = ve;
    options.colorfx = '';
    options.start_points = start_points;
    col = Q(:,1);
    if strcmp(test, 'voronoi')
        col(1) = nmax;
    end
    plot_fast_marching_mesh(vertex,faces, col, [], options);
    colormap(cm);
    saveas(gcf, [rep name '-' test '-' num2string_fixeddigit(i,2) '.png'], 'png');
    % update positions
    if strcmp(test, 'lloyd')
        options.edges_id = edges_id;
        options.lambda = lambda;
        options.Q = Q;
        start_points = perform_lloyd_mesh(vertex,faces, start_points, options);
    else
        % add new points
        start_points = [start_points; round(rand(nadd,1)*(n-1)+1)];
    end
end
