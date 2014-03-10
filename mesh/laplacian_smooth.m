function [U,Uall] = laplacian_smooth(V,F,L_method,b,lambda,method,S,max_iter)
  % LAPLACIAN_SMOOTH smooth a mesh using implicit/explicit laplacian smoothing
  %
  % [U] = laplacian_smooth(V,F)
  % [U] = laplacian_smooth(V,F,L_method,b,lambda,method,S)
  % [U,Uall] = laplacian_smooth(V,F,L_method,b,lambda,method,S,max_iter)
  % 
  % Inputs:
  %   V  #V x 3 matrix of vertex coordinates
  %   F  #F x 3  matrix of indices of triangle corners
  %   L_method  method for laplacian
  %      'uniform'
  %      'cotan'
  %   b  list of indices of fixed vertices
  %   lambda  diffusion speed parameter {0.1}
  %   method  method to use:
  %     'implicit' (default)
  %     'explicit'
  %   S  scalar fields to smooth (default V)
  % Outputs:
  %   U  #V x 3 list of new vertex positions
  %   Uall  #V x 3 x iters list of new vertex positions for each iteration
  %   

  % number of vertices
  n = size(V,1);

  % nuymber of dimensions
  dim = size(V,2);


  if(~exist('L_method','var'))
    if is_planar(V)
      L_method = 'uniform';
    else
      L_method = 'cotan';
    end
  end

  if(~exist('lambda','var')) 
    lambda = 0.1; 
  end

  if(~exist('b','var'))
    b = [];
  end

  % bulid sparse identity matrix
  I = speye(n,n);

  if(~exist('method','var'))
    method = 'implicit';
  end


  h = avgedge(V,F);
  if(~exist('tol','var'))
    tol = 0.001;
  end

  if(~exist('max_iter','var'))
    max_iter = 1000;
  end

  if(~exist('S','var'))
    S = V;
  end


  % only compute uniform laplacain once
  if strcmp(L_method,'uniform')
    % planar meshes should use uniform laplacian
    A = adjacency_matrix(F);
    L = A - diag(sum(A));
  end

  % place for factorization and symmtery flag used by min_quad_with_fixed
  P = [];
  sym = [];

  iter = 0;
  U = S;
  U_prev = S;
  if nargout >= 2
    Uall = [];
  end 

  % recompute laplacian
  if strcmp(L_method,'cotan')
    % other 3D meshes should use cotangent laplacian
    L = cotmatrix_embedded(V,F);
    %error
  end

  while( iter < max_iter && (iter == 0 || max(abs(U(:)-U_prev(:)))>tol*h))
    U_prev = U;

    switch method
    case 'implicit'
      Q = (I-lambda*L);
      % could prefactor Q for 'uniform' case
      for d = 1:size(S,2)
        [U(:,d),P] = min_quad_with_fixed(Q*0.5,-U(:,d),b,S(b,d),[],[],P);
      end
    case 'explicit'
      Q = (I+lambda*L);
      U = Q * U;
      % enforce boundary
      U(b,:) = S(b,:);
    otherwise
      error(['' method ' is not a supported smoothing method']);
    end

    if nargout >= 2
      Uall = cat(3,Uall,U);
    end
    iter = iter + 1;
  end

  %[iter max_iter]
  %[max(abs(U(:)-U_prev(:))) tol*h]

end
