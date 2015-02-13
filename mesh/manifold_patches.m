function [C,A,uE2F,uE] = manifold_patches(F)
  % MANIFOLD_PATCHES Compute connected components of facets connected by
  % manifold edges.
  %
  % [C,A,uE2F,uE] = manifold_patches(F)
  %
  % Known bugs: This will detect a moebius strip as a single patch (manifold,
  % non-orientable) and also non-manfiold, yet orientable patches. 
  %
  % Q: Does this find exactly (manifold || orientable) patches?
  %
  % Inputs:
  %   F  #F by simplex-size list of facets
  % Outputs:
  %   C  #F list of component ids
  %   A  #F by #F adjacency matrix
  %

  % simplex size
  ss = size(F,2);
  assert(ss == 3);

  [A,uE2F,uE] = facet_adjacency_matrix(F,'ManifoldOnly',true);
  % Connected components are patches
  [~,C] = conncomp(A);

end
