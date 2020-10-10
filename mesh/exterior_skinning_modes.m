function [U,N,Y] = exterior_skinning_modes(V,D,H,M,k)
  % EXTERIOR_SKINNING_MODES Perform modal analysis on H,M within the subspace
  % where vertices in "linear blend skinning handle" selected by D of a mesh with
  % vertices V deform affinely.
  % 
  % [U,N,Y,H,M] = exterior_skinning_modes(V,D,H,M,k)
  %
  % Inputs:
  %   V  #V by dim
  %   D  #V by #H bone-handle incidence matrix D(i,j) = 1 means vertex i belongs
  %     to bone j
  %   H  #V*dim by #V*dim stiffness matrix
  %   M  #V*dim by #V*dim mass matrix
  %   k  number of eigen modes to compute
  % Outputs:
  %   U  #V*3 by k list of eigen modes in the "maximal" space
  %   N  #V*3 by #N subspace reduction matrix
  %   Y  #N by k list of subspace modes so that U = N * Y
  %
  % Example:
  %   H = arap_hessian(V,F);
  %   M = repdiag(massmatrix(V,F),3);
  %   [~,N,Y] = exterior_skinning_modes(V,D,H,M,20);
  %

  assert(max( sum(D,2) == 1));
  assert(all( any(D,1) ));
  n = size(V,1);
  dim = size(V,2);
  assert(size(H,1) == n*dim);
  assert(size(M,1) == n*dim);
  assert(k <= n);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Build the null-space so that V = N * Y
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  int = find(~any(D,2));
  b = find(any(D,2));
  ni = numel(int);
  m = size(D,2);
  S = sparse(lbs_matrix(V(b,:),D(b,:)));
  NI = speye(dim*ni,dim*ni+size(S,2));
  vec = @(X) X(:);
  b3 = vec(b+(0:dim-1)*n);
  int3 = vec(int+(0:dim-1)*n);
  N = sparse(int3,1:numel(int)*dim,1,n*dim,dim*ni+size(S,2));
  N(b3,:) = [sparse(dim*(n-ni),dim*ni) S];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Test that applying identity transformations to handles produces V
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Astack = repmat([eye(3,3) zeros(3,1)],[1 1 m]);
  %% collect transformations into column
  %A = reshape(permute(Astack,[3 1 2]),m*3*(3+1),1);
  %Y = [reshape(V(int,:),[],1);A];
  %tsurf(F,reshape(N*Y,[],3),'CData',1*any(D,2),fphong)
  %axis equal;
  %view(90,0);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Handle constrained Eigen Decomposition
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % wow, this symmetrization is important...
  smash = @(A,P) 0.5*(P'*A*P + (P'*A*P)');
  [Y,YD] = eigs(-0.5*smash(H,N),0.5*smash(M,N),k,'sm');

  U = N*Y;
end
