function [C] = adjacency_edge_cost_matrix(V,E)
  % ADJACENCY_EDGE_COST_MATRIX Build sparse  adjacency matrix from edge list or
  % face list
  % 
  % [C] = adjacency_edge_cost_matrix(V,E)
  % [C] = adjacency_edge_cost_matrix(V,F)
  % [C] = adjacency_edge_cost_matrix(V,T)
  %
  % Inputs:
  %  V  #V list of vertex positions
  %  E  #E by 2 edges list
  %  or 
  %  F  #F by 3 triangle list
  %  or 
  %  T  #F by 4 tet list
  % Outputs:
  %   C  #V by #V adjacency matrix weighted by edge norms
  %

  if(size(E,2)>2)
    F = E;
    E = edges(F);
  end

  % compute edge norms
  edge_norms = sqrt(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2));
  % number of vertices
  n = size(V,1);
  % build sparse adjacency matrix with non-zero entries indicated edge costs
  C = sparse([E(:,1);E(:,2)],[E(:,2);E(:,1)],repmat(edge_norms,[2 1]),n,n);
end
