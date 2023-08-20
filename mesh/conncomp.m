function [S,C] = conncomp(G)
  % CONNCOMP Drop in replacement for graphconncomp.m from the bioinformatics
  % toobox. G is an n by n adjacency matrix, then this identifies the S
  % connected components C. This is also an order of magnitude faster.
  %
  % [S,C] = conncomp(G)
  %
  % Inputs:
  %   G  n by n adjacency matrix, G(i,j) = G(j,i) ≠ 0 implies ij are connected.
  % Outputs:
  %   S  scalar number of connected components
  %   C  

  % Transpose to match graphconncomp
  G = G';

  A = G+speye(size(G));
  [p,~,r] = dmperm(A);
  S = numel(r)-1;
  C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
  C(p) = C;
end
