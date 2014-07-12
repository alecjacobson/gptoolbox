function [layers,interior,layershe] = getlayers(mesh, nlayers)
% GETLAYERS Get layers from a manifold mesh stored as a half-edge
% datastructure.
%
% [layers,interior] = getlayers(mesh, nlayers)
%
% Inputs:
%   mesh  mesh struct
%   nlayers  number of layers to get, >= 1
% Outputs:
%   layers  cell array of index sets for layers 1.. k, and index set for
%     interior layers are defined recursively as layer(1) = vertices on the
%     boundary layer(i) = vertices edge-adjacent to layer(i-1)
%   interior  list of indices to interior vertices
%   layershe  list of layers as half-edges
%
% See also: hds

    ind = 1:mesh.nhe;
    layers = {};
    % assume that no more than one layer boundary halfedge points to a vertex
    % layer(1)
    layers{1} = mesh.tip(mesh.bndhe);
    layershe{1} = mesh.bndhe;
    allprev_layers = layers{1};
    for i = 2:nlayers
        % layer{i} is the set difference between the set of all tail vertices
        % of halfedges with tips pointing to vertices in layer{i-1}
        % and all vertices in layers 0.. i-1
        layers{i} = setdiff(...
        mesh.tail(ind(ismember(mesh.tip,layers{i-1}))),...
        allprev_layers);
        % now find the halfedges of the layer;
        % these are defined as the boundary halfedges of the mesh we obtain
        % if we remove all faces with a vertex in layer{1}.. layer{i-1}
        % this is equivalent (I think) to the halfedges satisfying the following conditions:
        % 1) it is the halfedge of a triangle with exactly two vertices in layer{i}
        % 2) the  remaining vertex of the triangle is in layer{i-1},
        % 3) its tip and tail in layer{i}
        % 4) the triangle on the opposite side is in not in layer{i-1}

        % triangles with 2 vertices in layer{i} and a vertex outside
        tri_with_i_edge = ismember(mesh.F(1,:), layers{i}) + ...
                          ismember(mesh.F(2,:), layers{i}) + ...
                          ismember(mesh.F(3,:), layers{i}) == 2;
        tri_with_outside_vert = tri_with_i_edge & ...
                               (ismember(mesh.F(1,:), layers{i-1}) + ...
                                ismember(mesh.F(2,:), layers{i-1}) + ...
                                ismember(mesh.F(3,:), layers{i-1}) == 1);
        % now get the indices of halfedges on the boundary
        tri_ind = find(tri_with_outside_vert);
        % indices of all triangle halfedges
        he_ind = [tri_ind tri_ind+mesh.nfaces tri_ind+2*mesh.nfaces];
        % the halfedges we need have both the tip and tail in layer{i}
        he_ind_layer_i  = he_ind(ismember(mesh.tip(he_ind),layers{i}) & ...
                             ismember(mesh.tail(he_ind),layers{i}) );
        % he pointing to the opposite vertex of the triangle on the other side
        he_ind_opp_vertex = mesh.next(mesh.opp(he_ind_layer_i));
        layershe{i} = he_ind_layer_i(~ismember(mesh.tip(he_ind_opp_vertex),layers{i-1}));
        allprev_layers = union(layers{i},allprev_layers);
    end
    interior = setdiff(1:mesh.nvert,allprev_layers);
end
