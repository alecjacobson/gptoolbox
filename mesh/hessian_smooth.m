function [U,Usteps] = hessian_smooth(V,F,varargin)
  % HESSIAN_SMOOTH smooth a planar mesh using implicit/explicit hessian
  % smoothing by minimizing the integrated squared hessian energy
  % see Stein et al. 2018, http://www.cs.columbia.edu/cg/hessians/
  %
  % [U] = hessian_smooth(V,F)
  % [U] = hessian_smooth(V,F,b,lambda,method,S)
  % [U,Usteps] = hessian_smooth(V,F,'b',b,'Lambda','S',f)
  % 
  % Inputs:
  %   V  #V x 2 matrix of vertex coordinates
  %   F  #F x 3  matrix of indices of triangle corners
  % Optional:
  %   b  list of indices of fixed vertices
  %   Lambda  diffusion speed parameter {0.1}
  %   Method  method to use:
  %     'implicit' (default)
  %     'explicit'
  %     'limit'
  %   S  scalar fields to smooth (default V)
  %   MaxIter  maximum number of iterations to solve
  %   MaxDiff  minimum difference between consecutive iterations for
  %            convergence
  % Outputs:
  %   U  smoothed function values
  %   Usteps  list of smoothed function values for each iteration
  %   
  % See also: hessian_squared, laplacian_smooth
  %
  
  % defaults
  lambda = 0.1;
  method = 'implicit';
  b = [];
  S = V;
  max_iter = 1000;
  max_diff = 1e-13;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'b','Lambda','Method','S','MaxIter','MaxDiff'}, ...
    {'b','lambda','method','S','max_iter','max_diff'});
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
  H = hessian_squared(V,F);
  M = massmatrix(V,F,'barycentric');
  I = speye(size(H));
  
  if strcmp(method, 'limit')
    for d = 1:size(S,2)
      U(:,d) = min_quad_with_fixed(0.5*H,[],b,S(b,d));
    end
    return
  end
  
  if nargout > 1
    Usteps = zeros([size(S) max_iter]);
  end
  U = S;
  iter = 1;
  while true
    if nargout > 1
      Usteps(:,:,iter) = U;
    end
    
    U_prev = U;
    switch method
    case 'implicit'
      Q = (M-lambda*H);
      for d = 1:size(S,2)
        U(:,d) = min_quad_with_fixed(Q*0.5,-U(:,d),b,S(b,d),[],[]);
      end
    case 'explicit'
      Q = (I+lambda*H);
      U = Q * U;
      % enforce boundary
      U(b,:) = S(b,:);
    otherwise
      error(['' method ' is not a supported smoothing method']);
    end
    
    % Use difference from previous as stopping criterion
    d = trace(((U-U_prev)'*M*(U-U_prev)).^2);
    if d < max_diff
      warning('converged...');
      break;
    end
    if iter >= max_iter
      warning('Max iterations (%d) exceeded without convergence',max_iter);
      break;
    end
    iter = iter + 1;
  end

  if nargout > 1
    Usteps = Usteps(:,:,1:iter);
  end

end

