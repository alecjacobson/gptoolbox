function [layershe] = ...
  getlayershe(mesh, outside_region_of_interest, nlayers)
% GETLAYERSHE Does the same as getlayers except instead of finding layers by
% shrinking from boundary it expects mesh.layers to already contain the proper
% indices. So that this function just sets layershe properly
% 
% [layers,interior] = getlayers(mesh, nlayers)
%
% Inputs:
%   mesh  mesh half-edge struct
%   nlayers  number of layers to get, >= 1
% Outputs:
%   lyaershe  cell array of index sets for layers 1.. k, and index set for
%     interior layers are defined recursively as layer(1) = vertices on the
%     boundary layer(i) = vertices edge-adjacent to layer(i-1)

    ind = 1:mesh.nhe;
    for i = 1:nlayers
        % layer{i} is the set difference between the set of all tail vertices
        % of halfedges with tips pointing to vertices in layer{i-1}
        % and all vertices in layers 0.. i-1
        %layers{i} = setdiff(... 
        %mesh.tail(ind(ismember(mesh.tip,layers{i-1}))),...
        %allprev_layers);

        % now find the halfedges of the layer;
        % these are defined as the boundary halfedges of the mesh we obtain
        % if we remove all faces with a vertex in layer{1}.. layer{i-1}
        % this is equivalent (I think) to the halfedges satisfying the following conditions:
        % 1) it is the halfedge of a triangle with exactly two vertices in layer{i} 
        % 2) the  remaining vertex of the triangle is in layer{i-1},  
        % 3) its tip and tail in layer{i}
        % 4) the triangle on the opposite side is in not in layer{i-1}
        
        % triangles with 2 vertices in layer{i} and a vertex outside
        tri_with_i_edge = ismember(mesh.F(1,:), mesh.layers{i}) + ...
                          ismember(mesh.F(2,:), mesh.layers{i}) + ...
                          ismember(mesh.F(3,:), mesh.layers{i}) == 2;


        if(i==1)
          % TODO this is not catching case when layers{1} has true boundary hes
          tri_with_outside_vert = tri_with_i_edge & ...
                                (ismember(mesh.F(1,:),  outside_region_of_interest) + ...
                                  ismember(mesh.F(2,:), outside_region_of_interest) + ...
                                  ismember(mesh.F(3,:), outside_region_of_interest) == 1);
        else
          tri_with_outside_vert = tri_with_i_edge & ...
                                (ismember(mesh.F(1,:),  mesh.layers{i-1}) + ...
                                  ismember(mesh.F(2,:), mesh.layers{i-1}) + ...
                                  ismember(mesh.F(3,:), mesh.layers{i-1}) == 1);
        end

        % now get the indices of halfedges on the boundary
        tri_ind = find(tri_with_outside_vert);
        % indices of all triangle halfedges
        he_ind = [tri_ind tri_ind+mesh.nfaces tri_ind+2*mesh.nfaces];
        % the halfedges we need have both the tip and tail in layer{i}
        he_ind_layer_i  = he_ind(ismember(mesh.tip(he_ind),mesh.layers{i}) & ...
                             ismember(mesh.tail(he_ind),mesh.layers{i}) );
        % he pointing to the opposite vertex of the triangle on the other side 
        he_ind_opp_vertex = mesh.next(mesh.opp(he_ind_layer_i));
        if(i==1)
          % TODO this is not catching case when layers{1} has true boundary hes
          layershe{i} = he_ind_layer_i(~ismember(mesh.tip(he_ind_opp_vertex),outside_region_of_interest));
        else
          layershe{i} = he_ind_layer_i(~ismember(mesh.tip(he_ind_opp_vertex),mesh.layers{i-1}));
        end
    end
end
    
    
    
       
       
       
       
