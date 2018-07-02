function [VV] = vv(V,F)
  warning('vv is deprecated, use adjacency_matrix instead.');
  VV = adjacency_matrix(F);
end
