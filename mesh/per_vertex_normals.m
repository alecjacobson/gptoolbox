function N = per_vertex_normals(V,F)
  % PER_VERTEX_NORMALS  Compute per-vertex normals over a mesh (V,F)
  %
  % N = per_vertex_normals(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  % Outputs:
  %   N  #V by 3 list of vertex normals
  %

  dblA = doublearea(V,F);
  A= sparse(F(:),repmat(1:size(F,1),1,3),repmat(dblA,3,1),size(V,1),size(F,1));
  FN = normalizerow(normals(V,F)+eps);
  N = normalizerow(A*FN);

end
