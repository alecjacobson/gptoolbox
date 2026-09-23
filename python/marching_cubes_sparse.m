function [V,F] = marching_cubes_sparse(S,GV,GI,isovalue)
  % [V,F] = marching_cubes(X,Y,Z,S)
  if nargin<4
    isovalue = 0;
  end

  assert(size(GV,2) == 3);
  assert(size(GI,2) == 8);
  GI_minus_1 = GI - 1; % Convert to zero-based indexing
  Svec = reshape(permute(S,[2 1 3]),[],1);
  res = pyrunfile('marching_cubes_sparse.py', 'res', S=Svec, GV=GV,GI=GI_minus_1, isovalue=isovalue);
  V = double(res.V);
  F = double(res.F)+1;
end

