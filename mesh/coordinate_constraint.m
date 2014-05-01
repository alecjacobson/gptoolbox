function [Alp,blp] = coordinate_constraint(V,C)
  % COORDINATE_CONSTRAINT  Build a "coordinate constraint" matrix for enforcing
  % linear precision in generalized barycentric coordinates or skinning
  % weights.
  %
  % [Alp,blp] = coordinate_constraint(V,C)
  % 
  % Inputs:
  %   V  #V by dim list of sample locations
  %   C  #C by dim list of cage/control vertex positions
  % Outputs:
  %   Alp  #V*dim by #V*c sparse matrix of constraint coefficients
  %   blp  #V*dim list of right hand sides
  %
  % Example:
  %   % fine mesh in (V,F) cage in (C,CF), boundary conditions in (b,bc)
  %   [Alp,blp] = coordinate_constraint(V,C)
  %   W = min_quad_with_fixed(Q,[],b,bc,Alp,blp)
  %   max(abs(Alp*W(:) - blp))
  % 

  % dimensions
  n = size(V,1);
  d = size(C,2);
  assert(d == size(V,2));
  c = size(C,1);
  AlpI = reshape(repmat(1:n,c,1),n*c,1);
  AlpJ = reshape(bsxfun(@plus,repmat((1:c)-1,n,1)*n,(1:n)')',n*c,1);
  %% Handle each dimension
  AlpI = reshape(bsxfun(@plus,AlpI,((1:d)-1)*n),n*c*d,1);
  AlpJ = repmat(AlpJ(:),d,1);
  AlpV = reshape(repmat(C,n,1),n*c*d,1);
  Alp = sparse(AlpI(:),AlpJ(:),AlpV(:),n*d,c*n);
  blp = V(:);

end
