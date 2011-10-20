function [in,on,B,L] = in_mesh(V,F,Q)
  % IN_MESH test whether a list of points are in a given mesh
  %  
  % [in] = in_mesh(V,F,Q)
  % [in,on,B,L] = in_mesh(V,F,Q)
  % 
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of face indices
  %   Q  #Q by dim list of query points
  % Outputs:
  %   in #Q list of flags revealing whether queries are in (V,F)
  %   on #Q list of flags revealing whether queries are on boundary of (V,F)
  %   B  #B by 1 list of mesh outline edges 
  %   L  #loops+1 by 1 list of boundary loop start indices into B, the last
  %     entries is (by tradition) always the numel of B + 1
  %
  % See inpolygon
  %

  dim = size(V,2);
  % only works in 2D
  assert(dim == 2);

  % get ordered mesh outline
  [B,L] = ordered_outline(F);

  % don't know how to handle holes / multiple components
  assert(numel(L) == 2);

  % only have one loop
  l = 1;
  [in, on] = inpolygon( ...
    Q(:,1), ...
    Q(:,2), ...
    V([B(L(l):(L(l+1)-1)) B(L(l))],1), ...
    V([B(L(l):(L(l+1)-1)) B(L(l))],2));

end
