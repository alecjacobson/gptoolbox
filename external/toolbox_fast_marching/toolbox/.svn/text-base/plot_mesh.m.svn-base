function h = plot_mesh(vertex,face,options)

% plot_mesh - plot a 3D mesh.
%
%   plot_mesh(vertex,face, options);
%
%   'options' is a structure that may contains:
%       - 'normal' : a (nvertx x 3) array specifying the normals at each vertex.
%       - 'edge_color' : a float specifying the color of the edges.
%       - 'face_color' : a float specifying the color of the faces.
%       - 'face_vertex_color' : a color per vertex or face.
%       - 'vertex'
%
%   See also: mesh_previewer.
%
%   Copyright (c) 2004 Gabriel Peyré


if nargin<2
    error('Not enough arguments.');
end
if nargin<3
    options.null = 0;
end


% can flip to accept data in correct ordering
if (size(vertex,1)==3 || size(vertex,1)==2) && size(vertex,2)~=3
    vertex = vertex';
end
if size(face,1)==3 && size(face,2)~=3
    face = face';
end

if size(face,2)~=3 || (size(vertex,2)~=3 && size(vertex,2)~=2)
    error('face or vertex does not have correct format.');
end


if ~isfield(options, 'normal')
    options.normal = [];
end
normal = options.normal;

if ~isfield(options, 'face_color')
    options.face_color = 0.7;
end
face_color = options.face_color;

if ~isfield(options, 'edge_color')
    options.edge_color = 0;
end
edge_color = options.edge_color;

if ~isfield(options, 'face_vertex_color')
    options.face_vertex_color = zeros(size(vertex,1),1);
end
face_vertex_color = options.face_vertex_color;


if isempty(face_vertex_color)
    h = patch('vertices',vertex,'faces',face,'facecolor',[face_color face_color face_color],'edgecolor',[edge_color edge_color edge_color]);
else
    nverts = size(vertex,1);
    % vertex_color = rand(nverts,1);
    if size(face_vertex_color,1)==size(vertex,1)
        shading_type = 'interp';
    else
        shading_type = 'flat';
    end
    h = patch('vertices',vertex,'faces',face,'FaceVertexCData',face_vertex_color, 'FaceColor',shading_type);
end
colormap gray(256);
lighting phong;
camlight infinite; 
camproj('perspective');
axis square; 
axis off;

if ~isempty(normal)
    % plot the normals
    hold on;
    quiver3(vertex(:,1),vertex(:,2),vertex(:,3),normal(:,1),normal(:,2),normal(:,3),0.8);
    hold off;
end

axis tight;
axis equal;
cameramenu;