function [E,C] = nonmanifold_edges(F)
  % NONMANIFOLD_EDGES List of non-manifold edges
  %
  % Inputs:
  %   F  #F by dim=3 list of facet indices
  % Outputs:
  %   E  #E by 2 list of nonmanifold edges
  %   C  #E by 1 list of unsigned Counts (>2)
  %   %pC  #E by 1 list of positive Counts
  %   %nC  #E by 1 list of negative Counts
  %

  allE = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
  sortallE = sort(allE,2);
  sC = sparse(sortallE(:,1),sortallE(:,2),(allE(:,1)<allE(:,2))*2-1);
  C = sparse(sortallE(:,1),sortallE(:,2),1);
%  assert(nnz(tril(C)) == 0);
% assert(nnz(tril(sC)) == 0);
  [EI,EJ,C] = find((C>2 | abs(sC)>1).* C);
  E = [EI EJ];
end
