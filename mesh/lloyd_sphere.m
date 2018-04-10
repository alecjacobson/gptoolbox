function [N,F] = lloyd_sphere(n,varargin)
  % LLOYD_SPHERE Construct a triangle mesh of the sphere with n reasonably well
  % distributed vertices.
  %
  % Inputs:
  %   n  number of points
  % Outputs:
  %   N  #N by 3 list of vertex positions
  %   F  #F by 3 list of face indices into N
  %

  %N = randsphere(n,'Method','trig');
  %N = normalizerow(N.^3);

  % default values
  subdivision_method = 'upsample';
  max_iters = 100000;
  N = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'InitialGuess','MaxIters','SubdivisionMethod'}, ...
    {'N','max_iters','subdivision_method'});
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

  switch subdivision_method
  case 'sqrt3'
    ssn = ceil(-(log(3)-log(3/10*n-3/5))/log(3));
  case {'upsample','loop'}
    ssn = ceil(log2(1/10*(10*(n)-20)^(1/2)));
  otherwise
    error('Unsupported subdivision_method %s',subdivision_method);
  end
  if isempty(N)
    N = subdivided_sphere(ssn, 'SubdivisionMethod',subdivision_method);
  end

  N = N(1:n,:);
  F = convhulln(N);
  if n <= 12
    % nothing more can be done
    return;
  end

  M = massmatrix(N,F,'voronoi');
  vis = false;
  if vis
    t = tsurf(F,N);
    caxis auto;
    set(t,'CData',full(diag(M)));
    colorbar;
    axis equal;
  end
  for iter = 1:max_iters
    A = adjacency_matrix(F);
    A = A*M;
    %A = bsxfun(@rdivide,A,sum(A,2));
    A = spdiags (1./sum (A,2), 0, size(A,1), size(A,1)) * A ;
    N_prev = N;
    N = A*N;
    % subtract off center of mass  (needed for small n)
    N = bsxfun(@minus,N,diag(M)'*N./sum(diag(M)));
    N = normalizerow(N);
    F = fliplr(convhulln(N));
    M = massmatrix(N,F,'voronoi');
    er = trace((N-N_prev)'*M*(N-N_prev));
    if vis
    set(t,'Vertices',N,'Faces',F, ...
      ... 'CData',full(diag(M)));
      'CData',full(sum(adjacency_matrix(F),2)),'FaceLighting','phong','FaceColor','interp');
    axis equal;
    drawnow;
    title(sprintf('%g',er));
    end
    if er < 1e-07
      break;
    end
  end
end
