function U = conformalized_mean_curvature_flow(V,F,varargin)
  % CONFORMALIZED_MEAN_CURVATURE_FLOW Flow a surface according to "Can mean
  % curvature flow be made non-singular?" [Kazhdan et al. 2012]
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   Optional:
  %     'MaxIter' followed by maximum number of iterations {100}
  %     'delta' followed by delta value, should roughly be in range
  %       [1e-13,1e13] {1}
  %     'LaplacianType' followed by 'cotangent' of 'uniform'.
  % Outputs:
  %   U  #V by dim list of new vertex positions
  %

  % default values
  delta = 1;
  max_iter = 100;
  laplacian_type = 'cotangent';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MaxIter','delta','LaplacianType'}, ...
    {'max_iter','delta','laplacian_type'});
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

  switch laplacian_type
  case 'cotangent'
    L = cotmatrix(V,F);
  case 'uniform'
    A = adjacency_matrix(F);
    L = A - diag(sparse(sum(A,2)));
  end

  U = V;
  iter = 1;
  while true
    U_prev = U;
    % 'full' seems slight more stable than 'barycentric' which is more stable
    % than 'voronoi'
    M = massmatrix(U,F,'barycentric');
    U = (M-delta*L)\(M*U);
    area = sum(doublearea(U,F)*0.5);
    c = sum(bsxfun(@times,0.5*doublearea(U,F)/area,barycenter(U,F)));
    U = bsxfun(@minus,U,c);
    U = U/sqrt(area);
    % Use difference from previous as stopping criterion
    % Better would be to look for convergence while factoring out Möbius
    % transformation.
    % Q: Stop when no change in angles?
    d = trace(((U-U_prev)'*M*(U-U_prev)).^2);
    if d < 1e-13
      break;
    end
    % Volume of unit area sphere: pi^-0.5/6
    [c,vol] = centroid(U,F);
    %tsurf(F,U);
    %%[vol pi^-0.5/6]
    %title(sprintf('%g',sum(doublearea(U,F)*0.5)));
    %drawnow;
    if iter > max_iter
      break;
    end
    iter = iter + 1;
  end
end
