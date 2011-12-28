function S_V = moveEV(V, E, S_E)

A = adjacency_matrix(E);
num_neighbors = sum(A,2);

S_V = zeros(size(V,1),1);
S_V(E(:,1)) = S_V(E(:,1)) + S_E;
S_V(E(:,2)) = S_V(E(:,2)) + S_E;

S_V = 0.5* S_V ./ num_neighbors;

end