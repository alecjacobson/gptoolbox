function [Q,H,M4,D,A,G,int] = hessian_squared(V,F,varargin)
  % HESSIAN_SQUARED Construct a matrix to compute the integrated squared
  % Hessian of a function over a triangle mesh using the mixed finite element
  % method.
  % 
  % [Q] = hessian_squared(V,F)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle mesh indices
  %   extraint #any boundary vertices that are explicitly to be treated as
  %   the interior
  % Outputs:
  %   Q  #V by #V sparse matrix so that X'*Q*X measures the integrated squared
  %     Hessian energy of a scalar function X
  %
  
  extraints = [];
  if(length(varargin)>0)
      extraints = varargin{1};
  end

  % Number of faces
  m = size(F,1);
  % Number of vertices
  n = size(V,1);
  % Number of dimension
  dim = size(V,2);
  m = size(F,1);
  n = size(V,1);
  G = grad(V,F);
  assert(size(G,1) == dim*m,'Gradient should equal dim*m');
  M = massmatrix(V,F);
  M4 = repdiag(M,dim^2);
  % Block transpose of G
  GG = sparse(m,0);
  for d = 1:dim
    GG = [GG G((d-1)*m+(1:m),:)];
  end
  D = repdiag(GG,dim);

  %% This is not necessary. The gradient operator G is guaranteed by
  %%construction to generate vectors in the face-plane. Therefore, even if the
  %%matrix divergence lies off this plane, the dot product with the gradient will
  %%precisely ignore the off-plane component.
  % project_div = true;
  % if project_div
  %   warning('projectin''');
  %   N = normalizerow(normals(V,F));
  %   DN = sparse(repmat(1:m,3,1)',reshape(1:m*3,m,3),N,m,m*3);
  %   D = (speye(m*3)-DN'*DN)*D;
  % end


  switch size(F,2)
  case 3
    A = repdiag(diag(sparse(doublearea(V,F)*0.5)),dim);
  case 4
    A = repdiag(diag(sparse(volume(V,F))),dim);
  end
  H = D'*A*G;

  b = unique(boundary_faces(F));
  int = setdiff(1:n,b);
  int = [int(:); extraints];
  int4 = (0:(dim^2-1))*n + int;
  notint4 = setdiff(1:(n*dim^2), int4);
  %M4 = M4(int4,int4);
  %H = H(int4,:);
  H(notint4,:) = 0;
  Q = H'*(M4\H);
end
