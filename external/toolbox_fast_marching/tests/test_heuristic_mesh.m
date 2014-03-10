% test for heuristically driven computations on meshes

n = 300;
repimg = 'results/heuristic-mesh/';
if exist(repimg)~=7
    mkdir(repimg);
end


path(path, '../toolbox_graph/off');
rep = '../toolbox_graph/off/';
name = 'nefertiti.off';
name = 'beetle.off';
name = 'fandisk.off';
name = 'david50kf.off';
name = 'skull.off';
[vertex,faces] = read_mesh([rep name]);

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end
if size(faces,1)>size(faces,2)
    faces = faces';
end

save_images = 1;

nverts = max(size(vertex));
start_points = [1];
end_points = [2000];
start_points = start_points(:);

% compute distance to starting point
[H,Z,Q] = perform_fast_marching_mesh(vertex, faces, end_points);


options.end_points = end_points(:);

lambda = 1;
options.heuristic = lambda*H(:);
[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, start_points, options);

col = D;
col(col==Inf) = 0;
options.face_vertex_color = col;
clf;
hold on;
plot_mesh(vertex, faces, options);
h = plot3(vertex(1,start_points),vertex(2,start_points), vertex(3,start_points), 'r.');
set(h, 'MarkerSize', 25);
h = plot3(vertex(1,end_points),vertex(2,end_points), vertex(3,end_points), 'b.');
set(h, 'MarkerSize', 25);
hold off;
colormap jet(256);
view(0,0);
shading interp;
% lighting none;
camlight;