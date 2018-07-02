function [Fp, Fi] = tt(F)
  warning('tt is deprecated, use triangle_triangle_adjacency instead.');
  [Fp, Fi] = triangle_triangle_adjacency(F);
end
