% test for front propagation on 3D meshes
%
%   Copyright (c) 2007 Gabriel Peyre

path(path, 'toolbox/');
path(path, '../toolbox_graph_data/');
path(path, '../toolbox_graph_data/off/');


rep = ['results/geodesic-mesh/'];
if not(exist(rep))
    mkdir(rep);
end


disp('Loading mesh.');
if not(exist('name'))
    name = 'bunny';
    name = 'elephant-50kv';
    name = 'david50kf';
    name = 'david_head';
    name = 'hand';
end
options.name = name;
[vertex,faces] = read_mesh(name);
clf;
plot_mesh(vertex, faces,options);
saveas(gcf, [rep name '-mesh.png'], 'png');

nverts = max(size(vertex));
if not(exist('nstart'))
    nstart = 1;
end
start_points = [1 round(nverts/2)];
start_points = floor(rand(nstart,1)*nverts)+1;
start_points = start_points(:);
options.end_points = [];
% enforce a nice starting position for the first point
switch name
    case 'bunny'
        start_points(1) = 8900;
    case 'elephant-50kv'
        start_points(1) = 24575;
    case 'david50kf'
        start_points(1) = 20361;
    case 'david_head'
        start_points(1) = 18080;
    case 'hand'
        start_points(1) = 29719;
end

disp('Performing propagation.');
[D,S,Q] = perform_fast_marching_mesh(vertex, faces, start_points, options);


% compute geodesics
npaths = max(20,4*nstart);
% npaths = 1;
[tmp,I] = sort( D(:) ); I = I(end:-1:1); I = I(1:round(nverts*1));
end_points = floor( rand(npaths,1)*(length(I)-1) )+1;
end_points = I(end_points);
% [tmp,I] = sort( D(:) ); end_points(1) = I(end);

options.v2v = compute_vertex_ring(faces);
options.e2f = compute_edge_face_ring(faces);

disp('Extracting geodesics');
options.method = 'discrete';
options.method = 'continuous';
paths = compute_geodesic_mesh(D,vertex,faces, end_points, options);
    
options.colorfx = 'equalize';
plot_fast_marching_mesh(vertex,faces, D, paths, options);
saveas(gcf, [rep name '-geodesics-' num2str(nstart) '.png'], 'png');
