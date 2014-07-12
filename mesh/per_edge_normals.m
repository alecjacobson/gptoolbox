function [N,E,EMAP] = per_edge_normals(V,F)
  % PER_EDGE_NORMALS  Compute per-edge normals over a mesh (V,F)
  %
  % N = per_edge_normals(V,F)
  % [N,E,EMAP] = per_edge_normals(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  % Outputs:
  %   N  #E by 3 list of edge normals 
  %   E  #E by 2 list of edge indices
  %   EMAP  #F*3 by 1 unique map of edges
  %

  m = size(F,1);
  allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
  % Map duplicate edges to first instance
  [E,~,EMAP] = unique(sort(allE,2),'rows');
  % Map each face-edge to a unique edge
  F2E = repmat(1:m,3,1)';
  dblA = doublearea(V,F);
  % #E by #F*3
  A = sparse(EMAP,F2E(:),repmat(dblA,3,1),size(E,1),m);
  FN = normalizerow(normals(V,F)+eps);
  N = normalizerow(A*FN);
end
