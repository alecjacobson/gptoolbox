function [A] = adjacency_matrix(E)
  % ADJACENCY_MATRIX Build sparse adjacency matrix from edge list or face list
  % 
  % [A] = adjacency_matrix(E)
  % [A] = adjacency_matrix(F)
  % [A] = adjacency_matrix(T)
  %
  % Inputs:
  %   E  #E by 2 edges list
  %   or 
  %   F  #F by 3 triangle list
  %   or 
  %   T  #F by 4 tet list
  % Outputs:
  %   A  #V by #V adjacency matrix (#V = max(E(:)))
  %    
  % See also: facet_adjacency_matrix
  %

  if size(E,2)>2
    F = E;
    E = edges(F);
  end

  A = sparse([E(:,1) E(:,2)],[E(:,2) E(:,1)],1);
end
