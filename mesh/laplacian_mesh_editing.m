function U = laplacian_mesh_editing(V,F,b,bc)
  % LAPLACIAN_MESH_EDITING deform a mesh using laplacian
  % mesh editing scheme, "Laplacian Surface Editing" by Olga Sorkine,
  %
  % U = laplacian_mesh_editing(V,F,b,bc)
  %
  % Inputs:
  %   V  #V by dim list of rest domain positions
  %   F  #F by 3 list of triangle indices into V
  %   b  #b list of indices of constraint (boundary) vertices
  %   bc  #b by dim list of constraint positions for b
  % Outputs:
  %   U  #V by dim list of new positions

  % number of vertices
  n = size(V,1);
  assert(max(b) <= n);
  assert(min(b) >= 1);
  % dimension
  dim = size(V,2);
  assert(dim == size(bc,2));

  [A, rhs] = laplacian_editing_system(V,F,b,bc);
  U = zeros(n*dim,1);
  indices = 1:n;
  interior = indices(~ismember(indices,b));
  interior_dim = repmat(interior,1,dim);
  interior_dim = interior_dim + ...
    reshape(repmat(((1:dim)-1)*n,numel(interior),1),size(interior_dim));
  b_dim = repmat(b,1,dim);
  b_dim = b_dim + ...
    reshape(repmat(((1:dim)-1)*n,numel(b),1),size(b_dim));
  U(interior_dim) = A \ rhs;
  U(b_dim) = bc(:);
  U = reshape(U,size(V));
end
