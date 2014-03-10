% test for farthest point sampling on 3D meshes

n = 300;
repimg = 'results/farthest-sampling-mesh/';
if exist(repimg)~=7
    mkdir(repimg);
end


path(path, '../toolbox_graph_data/off/');
name = 'nefertiti';
name = 'beetle';
name = 'fandisk';
name = 'david50kf';
name = 'brain';
name = 'pelvis';
[vertex,faces] = read_mesh(name);
options.name = name;

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end
if size(faces,1)>size(faces,2)
    faces = faces';
end

if strcmp(name, 'brain')
    % remove half of the vertices
    [tmp,I] = sort(vertex(3,:), 2, 'descend');
    J = reverse_permutation(I);
    vertex = vertex(:,I);
    faces = J(faces);    
    n = size(vertex,2); n = round(n/2);
    vertex = vertex(:,1:n);
    faces = faces(:,sum(faces>n)==0);
end


save_images = 1;

% plot sampling location
i = 0;
landmark = [];
for nbr_landmarks = [200 500 1000 2000 5000 10000]  % 100:50:500
    i = i+1;
    
    disp('Perform farthest point sampling.');
    landmark = perform_farthest_point_sampling_mesh( vertex,faces, landmark, nbr_landmarks-length(landmark) );
    
    % compute the associated triangulation
    [D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmark);
    [vertex_voronoi,faces_voronoi] = compute_voronoi_triangulation_mesh(Q,vertex,faces);
    options.method = 'fast';
    options.method = 'slow';
    faces_voronoi = perform_faces_reorientation(vertex_voronoi,faces_voronoi, options);
    
    % display
    col = D; col(col==Inf) = 0;
%    options.face_vertex_color = col;
    clf;
    hold on;
    plot_mesh(vertex, faces, options);
    h = plot3(vertex_voronoi(1,:),vertex_voronoi(2,:),vertex_voronoi(3,:), 'r.');
    set(h, 'MarkerSize', 20);
    hold off;
    colormap gray(256);
    shading interp;   
    camlight;
    
    if save_images
        saveas(gcf, [repimg name '-sampling-' num2string_fixeddigit(nbr_landmarks,4) '.png'], 'png');
    end
           
    clf;
    plot_mesh(vertex_voronoi,faces_voronoi, options);
    shading faceted; lighting flat;
    camlight;
    
    if save_images
        saveas(gcf, [repimg name '-remeshing-' num2string_fixeddigit(nbr_landmarks,4) '.png'], 'png');
    end
end