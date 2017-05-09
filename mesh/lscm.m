function [U,Q] = lscm(V,F,b,bc,Aeq,Beq,varargin)
  % LSCM Compute Least Squares Conformal Mapping for mesh
  %
  % U = lscm(V,F,b,bc)
  % U = lscm(V,F,b,bc,Aeq,Beq)
  %
  % Inputs:
  %   V  #V by dim list of rest domain positions
  %   F  #F by 3 list of triangle indices into V
  %   b  #b list of indices of constraint (boundary) vertices
  %   bc  #b by 2 list of constraint positions for b
  %   Aeq   #Aeq by 2*#V matrix of linear equality constraints {[]}
  %   Beq   #Aeq vector of linear equality constraint right-hand sides {[]}
  %   Optional:
  %     'Method' followed by one of the following:
  %        'desbrun'  "Intrinsic Parameterizations of Surface Meshes" [Desbrun
  %          et al. 2002]
  %        'levy'  "Least Squares Conformal Maps for Automatic Texture Atlas
  %          Generation" [Lévy et al. 2002]
  %        {'mullen'}  "Spectral Conformal Parameterization" [Mullen et al. 2008]
  % Outputs:
  %   U  #V by 2 list of new positions
  %   Q  #V*2 by #*2  quadratic coefficients matrix
  %
  % Note: This is the same system as takeo_asap up to a factor of 2.5
  %
  % See also: arap, takeo_arap, takeo_asap
  %


  method = 'mullen';
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Method'}, ...
    {'method'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end


  if nargin<=4
    Aeq = [];
    Beq = [];
  end
  
  % number of vertices
  n = size(V,1);
  % number of triangles
  nt = size(F,1);
  % number of original dimensions
  dim = size(V,2);

  switch method
  case 'levy'
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
  
    % compute gradient matrix
    G = grad(V,F);
  
    % Extract each coordinate's block
    Gx = G(1:nt,:);
    Gy = G(nt+(1:nt),:);
  
    % Triangle areas
    TA = repdiag(diag(sparse(doublearea(V,F))/2),2);
  
    % Build quadratic coefficients matrix
    Q = [Gx -Gy;Gy Gx]'*TA*[Gx -Gy;Gy Gx];
  
    % solve
    U = min_quad_with_fixed(Q,zeros(2*n,1),[b b+n],bc(:));
    % reshape into columns
    U = reshape(U,n,2);
  case 'mullen'
    A = vector_area_matrix(F);
    L = repdiag(cotmatrix(V,F),2);
    Q = -L - 2*A;
    U = min_quad_with_fixed(Q,zeros(2*n,1),[b b+n],bc(:),Aeq,Beq);
    % reshape into columns
    U = reshape(U,n,2);
  case 'desbrun'
    error('not implemented.');
    % This implements the Dirichlet + Chi energies but usually when people
    % refer to the [Desbrun et al. 2002] paper they mean the [Mullen et
    % al. 2008] interpretation: zero chi energy + ""natural"" boundary
    % conditions.
    lambda = 1;
    mu = 1;
    L = cotmatrix(V,F);
    % chi energy
    C = cotangent(V,F);
    l = edge_lengths(V,F);
    X = sparse( ...
      F(:,[1 2 3 1 2 3]), ...
      F(:,[2 3 1 3 1 2]), ...
      C(:,[2 3 1 3 1 2])./ ...
      l(:,[3 1 2 2 3 1]), ...
      size(V,1),size(V,1));
    Q = repdiag(lambda*L+mu*X,2);
  otherwise
    error('unsupported method: %s',method);
  end
end
