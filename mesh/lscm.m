function [U] = lscm(V,F,b,bc)
  % LSCM Compute Least Squares Conformal Mapping for mesh
  %
  % U = lscm(V,F,b,bc)
  %
  % Inputs:
  %   V  #V by dim list of rest domain positions
  %   F  #F by 3 list of triangle indices into V
  %   b  #b list of indices of constraint (boundary) vertices
  %   bc  #b by 2 list of constraint positions for b
  % Outputs:
  %   U  #V by 2 list of new positions
  %
  % Note: This is the same system as takeo_asap up to a factor of 2.5
  %
  % See also: arap, takeo_arap, takeo_asap
  %
  
  % number of vertices
  n = size(V,1);
  % number of triangles
  nt = size(F,1);
  % number of original dimensions
  dim = size(V,2);

  % first need to convert each triangle to its orthonormal basis, if coming
  % from 3D
  assert(dim == 2);


  %% Indices of each triangle vertex, I, and its corresponding two neighbors, J
  %% and K
  %I = [F(:,1)];
  %J = [F(:,2)];
  %K = [F(:,3)];

  %X = [V(I,1) V(J,1) V(K,1)];
  %Y = [V(I,2) V(J,2) V(K,2)];

  %WRe = [X(:,3)-X(:,2) X(:,1)-X(:,3) X(:,2)-X(:,1)];
  %WIm = [Y(:,3)-Y(:,2) Y(:,1)-Y(:,3) Y(:,2)-Y(:,1)];

  %% sqrt root of twice the area of each triangle
  %dT = sqrt(doublearea(V,F));

  %% build M matrix, real and imaginary parts
  %II = [1:nt 1:nt 1:nt];
  %JJ = [I;J;K]';
  %VVRe = [WRe(:,1)./dT WRe(:,2)./dT WRe(:,3)./dT];
  %VVIm = [WIm(:,1)./dT WIm(:,2)./dT WIm(:,3)./dT];

  %WWRe = sparse(II,JJ,WRe,nt,n);
  %WWIm = sparse(II,JJ,WIm,nt,n);
  %% These look like blocks in the gradient matrix
  %MRe = sparse(II,JJ,VVRe,nt,n);
  %MIm = sparse(II,JJ,VVIm,nt,n);

  %% build A matrix
  %A = [MRe -MIm; MIm MRe];

  %% quadratic system matrix
  %Q = A'*A;

  % Or equivalently

  % We want that the Cauchy Riemann equations hold:
  %
  %  ∂u/∂x = ∂v/∂y and ∂u/∂y = -∂v/∂x
  %
  % Instead we satisfy this in a least-squares sense, minimizing
  %
  %  ∫ (∂u/∂x - ∂v/∂y)² + (∂u/∂y + ∂v/∂x)² dA
  %
  % On a mesh we know how to compute a [∂u/∂x;∂u/∂y] = G u = [Gx;Gy] u, so this
  % quadratic energy becomes U' [Gx -Gy;Gy Gx]'*TA*[Gx -Gy;Gy Gx] U = U' Q U
  % with TA being the triangle areas, stacking U = [u;v]
  % 
  % This skips the complex-number notation employed by "Least Squares Conformal
  % Maps for Automatic Texture Atlas Generation" [Lévy et al. 2002]
  % 

  % compute gradient matrix
  G = grad(V,F);

  % Extract each coordinate's block
  Gx = G(1:nt,:);
  Gy = G(nt+(1:nt),:);

  % Triangle areas
  TA = repdiag(diag(sparse(doublearea(V,F))/2),2);

  % Build quadratic coefficients matrix
  Q = [Gx -Gy;Gy Gx]'*TA*[Gx -Gy;Gy Gx];

  %% This is completely equivalent to [Mullen et al. 2008] "Spectral Conformal
  %% Parameterization"
  %Q = -(repdiag(cotmatrix(V,F),2) + 2*vector_area_matrix(F));

  % solve
  U = min_quad_with_fixed(Q,zeros(2*n,1),[b b+n],bc(:));
  % reshape into columns
  U = reshape(U,n,2);
end
