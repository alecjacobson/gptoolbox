function S_V = moveEV(V, E, S_E)
  % MOVEEV Seems to give each vertex the average value across its incident
  % edges?
  %
  % S_V = moveEV(V, E, S_E)
  %
  % Inputs:
  %   V  #V by dim mesh vertex positions?
  %   E  #E by 2 list of edge indices?
  %   S_E  #E by 1 scalar values on edges?
  % Outputs:
  %   S_V  #V by 1 scalar vavlues on vertices?

A = adjacency_matrix(E);
num_neighbors = sum(A,2);

S_V = zeros(size(V,1),1);
S_V(E(:,1)) = S_V(E(:,1)) + S_E;
S_V(E(:,2)) = S_V(E(:,2)) + S_E;

S_V = 0.5* S_V ./ num_neighbors;

end
