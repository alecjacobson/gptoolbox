function A = lischinski(L,varargin)
  % LISCHINSKI Construct the discrete laplacian for a given lightness image L
  % according to "Interactive Local Adjustment of Tonal Values" by [Lischinski
  % et al. 2006]
  %
  % A = lischinski(L)
  % A = lischinski(L,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   L  h by w lightness image
  %   Optional:
  %     'Lambda'  followed by lambda value {0.2}
  %     'Alpha'  followed by alpha value {1}
  %     'Epsilon'  followed by epsilon value {0.0001}
  % Outputs:
  %   A  h*w by h*w sparse laplacian
  %

  % default values
  lambda = 0.2;
  alpha = 1;
  epsilon = 0.0001;

  % Parse optional arguments
  v = 1;
  while(v<=numel(varargin))
    switch varargin{v}
    case 'Lambda'
      assert((v+1)<=numel(varargin));
      v = v+1;
      lambda = varargin{v};
    case 'Alpha'
      assert((v+1)<=numel(varargin));
      v = v+1;
      alpha = varargin{v};
    case 'Epsilon'
      assert((v+1)<=numel(varargin));
      v = v+1;
      epsilon = varargin{v};
    otherwise
      error(['Unsupported parameter: ' varargin{v}]);
    end
    v = v+1;
  end

  % width and height
  h = size(L,1);
  w = size(L,2);
  % number of pixels
  n = h * w;

  % indices of pixels
  indices = reshape(1:w*h,h,w);
  % down edges
  down_I = indices(1:end-1,1:end);
  down_J = indices(2:end,1:end);
  % across edges
  across_I = indices(1:end,1:end-1);
  across_J = indices(1:end,2:end);
  I = [down_I(:);across_I(:)];
  J = [down_J(:);across_J(:)];
  % build lower triangle
  A = sparse(I,J,-lambda*(abs(L(I)-L(J)).^alpha+epsilon).^-1,n,n);
  % matrix is symmetric
  A = A + A';
  % diagonal entires are -sum of off diagonals
  A = A - diag(sum(A,2));

  
end
