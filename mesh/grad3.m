function G = grad3(V,T)
  % GRAD3
  % 
  % G = grad(V,F)
  %
  % Compute the numerical gradient operator for tet meshes in 3d
  %
  % Inputs:
  %   V  #vertices by dim list of mesh vertex positions
  %   T  #elements by 3 list of mesh tet indices
  % Outputs:
  %   G  #elements*dim by #V Gradient operator
  %
  % Example:
  %   L = cotmatrix3(V,T);
  %   G  = grad3(V,T);
  %   vol = volume(V,T);
  %   GMG = -G'*repdiag(diag(sparse(vol)),3)*G;
  %

  % number of dimensions
  dim = size(V,2);
  assert(dim == 3);
  % number of vertices
  n = size(V,1);
  % number of elements
  m = size(T,1);
  % simplex size
  assert(size(T,2) == 4);

  % f(x) is piecewise-linear function:
  %
  % f(x) = ∑ φi(x) fi, f(x ∈ T) = φi(x) fi + φj(x) fj + φk(x) fk + φl(x) fl
  % ∇f(x) = ...                 = ∇φi(x) fi + ∇φj(x) fj + ∇φk(x) fk + ∇φl(x) fl
  %                             = ∇φi fi + ∇φj fj + ∇φk fk + ∇φl fl
  %
  % ∇φi = 1/hjk = Ajkl / 3V * (Facejkl)^perp
  %     = Ajkl / 3V * (Vj-Vk)x(Vl-Vk)
  %     = Ajkl / 3V * Njkl / ||Njkl||
  % 

  % get all faces
  F = [ ...
    T(:,1) T(:,2) T(:,3); ...
    T(:,1) T(:,3) T(:,4); ...
    T(:,1) T(:,4) T(:,2); ...
    T(:,2) T(:,4) T(:,3)];
  % compute areas of each face
  A = doublearea(V,F)/2;
  N = normalizerow(normals(V,F));

  % compute volume of each tet
  vol = volume(V,T);

  G = sparse( ...
    [0*m + repmat(1:m,1,4) ...
     1*m + repmat(1:m,1,4) ...
     2*m + repmat(1:m,1,4)], ...
    repmat([T(:,4);T(:,2);T(:,3);T(:,1)],3,1), ...
    repmat(A./(3*repmat(vol,4,1)),3,1).*N(:), ...
    3*m,n);

end
