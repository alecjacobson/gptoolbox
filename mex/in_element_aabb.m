% IN_ELEMENT_AABB Use an axis-aligned bounding box tree to determine if a set
% of points appears inside the elements of a mesh.
%
% I = in_element_aabb(V,Ele,Q);
% [I,bb_mins,bb_maxs,elements] = ...
%   in_element_aabb(V,Ele,Q,bb_mins,bb_maxs,elements);
%
% Inputs:
%   V  #V by dim list of mesh vertex positions. 
%   Ele  #Ele by dim+1 list of mesh indices into #V. 
%   Q  #Q by dim list of query point positions
%   Optional:
%     bb_mins  max_tree by dim list of bounding box min corner positions
%     bb_maxs  max_tree by dim list of bounding box max corner positions
%     elements  max_tree list of element or (not leaf id) indices into Ele
% Outputs:
%   I  #Q list of indices into Ele of first containing element (0 means no
%     containing element)
%   bb_mins (see optional input)
%   bb_maxs (see optional input)
%   elements (see optional input)
%
