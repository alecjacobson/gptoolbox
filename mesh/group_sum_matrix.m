function A = group_sum_matrix(G,k)
  % GROUP_SUM_MATRIX Builds a matrix A such that A*V computes the sum of
  % vertices in each group specified by G
  %
  % A = group_sum_matrix(G,k)
  % 
  % Inputs:
  %   G  #V list of group indices (1 to k) for each vertex, such that vertex i 
  %     is assigned to group G(i)
  %   k  #groups, default is max(G)
  % Outputs:
  %   A  #groups by #V sparse matrix such that A*V = group_sums
  %
  % See also: centroid_matrix


  % number of vertices
  n = numel(G);
  % number of groups
  if ~exist('k','var')
    k = max(G);
  end

  indices =  1:n;

  % builds A such that A(i,j) = 1 where i corresponds to group i and j
  % corresponds to vertex j
  A = sparse(G,indices,1,k,n);

end
