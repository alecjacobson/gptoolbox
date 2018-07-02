function [VT] = vertex_triangle_adjacency(F)
  % VERTEX_TRIANGLE_ADJACENCY Build a vertex-triangle adjacency matrix for a
  % triangle mesh.
  %
  % VT = vertex_triangle_adjacency(F)
  %
  % Input:
  %   F   #F x 3   list of face indices
  % Output:
  %   VT  #F x #V  sparse matrix, VT(i,j) is 1 iff vertex i is in face j

  i = (1:size(F,1))';
  j = F;
  VT = sparse([i i i],j,1);

end
