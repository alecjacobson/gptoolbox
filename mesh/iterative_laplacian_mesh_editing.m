function [U,preF] = iterative_laplacian_mesh_editing(varargin)
  % ITERATIVE_LAPLACIAN_MESH_EDITING deform a mesh using iterative laplacian
  % mesh editing scheme (described many places, for example "Handle-Aware
  % Isolines for Scalable Shape Editing")
  %
  % U = iterative_laplacian_mesh_editing(V,F,b,bc)
  %
  % Inputs:
  %   V  #V by dim list of rest domain positions
  %   F  #F by 3 list of triangle indices into V
  %   b  #b list of indices of constraint (boundary) vertices
  %   bc  #b by dim list of constraint positions for b
  %   Optional:
  %     'V0' #V by dim list of initial guess positions
  %       dim by dim by #C list of linear transformations initial guesses,
  %       optional (default is to use identity transformations)
  %     'Tol'
  %       stopping critera parameter. If variables (linear transformation matrix
  %       entries) change by less than 'tol' the optimization terminates,
  %       default is 0.75 (weak tolerance)
  %     'MaxIter'
  %       max number of local-global iterations, default is 10
  %     'PreF'
  %       prefactored systems for each dimension
  % Outputs:
  %   U  #V by dim list of new positions

  % parse input
  V = varargin{1};
  F = varargin{2};
  b = varargin{3};
  bc = varargin{4};

  % number of vertices
  n = size(V,1);
  assert(max(b) <= n);
  assert(min(b) >= 1);
  % dimension
  dim = size(V,2);
  assert(dim == size(bc,2));

  max_iterations = 100;
  tol = 0.001;
  preF = []; 

  ii = 5;
  while(ii <= nargin)
    switch varargin{ii}
    case 'V0'
      ii = ii + 1;
      assert(ii<=nargin);
      U = varargin{ii};
    case 'Tol'
      ii = ii + 1;
      assert(ii<=nargin);
      tol = varargin{ii};
    case 'MaxIter'
      ii = ii + 1;
      assert(ii<=nargin);
      max_iterations = varargin{ii};
    case 'PreF'
      ii = ii + 1;
      assert(ii<=nargin);
      if ~isempty(varargin{ii})
        preF = varargin{ii};
      end
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii + 1;
  end

  if(~exist('U','var') || isempty(U))
    U = V;
  end

  % Minimize the differences between the Laplacian coordinates before and after
  % editing in a least squares sense:
  % 
  %   arg min || LX - δ(X) ||²
  %        X
  %
  % Equivalent to solving the following in a least-squares sense:
  %
  % A * X = B(X)
  %
  % where A, laplace operator of original mesh; B(X) = δ(X), the laplacian
  % coodinates of the deformed mesh, that is, δ is the laplace operator of the
  % deformed mesh
  %
  % Minimize this by iterating two steps:
  %   Compute B0 = B(X0) = δ0(X0)
  %     that is build the laplace operator according to X0 and apply to X0
  %   Solve A * X = B0 in least squares sense
  %     that is:
  %     min (A * X - B0)' * (A * X - B0)
  %     min X' * A' * A * X - 2 * X' * A' * B0 + B0' * B0
  %     min X' * Q * X + X' * L + C
  %     with:
  %       Q = A' * A
  %       L = - 2 * A' * B0
  %     solve with:
  %       X = - Q \ ( L * X )
  %     
  %   All subject to position constraints on boundary (handles)
  %


  iteration = 0;
  U_prev = V;
  if ~isfield(preF,'mqwf')
    preF.LV = cotmatrix(V,F);
    preF.Q = preF.LV'*preF.LV;
    preF.mqwf = [];
    preF.h = avgedge(V,F);
  else
    preF.Q = [];
  end
  while( ...
    iteration < max_iterations && ...
    (iteration == 0 || max(abs(U(:)-U_prev(:)))>tol*preF.h))
    % remember last solution
    U_prev = U;

    % always enforce boundary conditions
    U(b,:)  = bc;

    % build laplace operator on deformed mesh
    LU = cotmatrix(U,F);
    % compute laplacian coordinates on deformed mesh
    % in 2D this is always zero
    B = LU * U;
   
    %E = sum(sum((preF.LV*U - LU * U).^2,2));

    % solve for new positions
    [U,preF.mqwf] = min_quad_with_fixed(preF.Q,-2*preF.LV'*B,b,bc,[],[],preF.mqwf);

    iteration = iteration + 1;
  end


end
