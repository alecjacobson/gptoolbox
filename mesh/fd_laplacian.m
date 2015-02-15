function L = fd_laplacian(varargin)
  % FD_LAPLACIAN  build a finite difference laplacian.
  %
  % L = fd_laplacian(h,w)
  % L = fd_laplacian([h,w])
  % 
  % Inputs:
  %   h  number of nodes along height 
  %   w  number of nodes along width
  % Outputs:
  %   L   h*w by h*w Laplacian with negative diagonal.
  % 
  function B = vec(A)
    B = A(:)';
  end

  switch nargin
  case 1
    h = varargin{1}(1);
    w = varargin{1}(2);
  case 2
    h = varargin{1};
    w = varargin{2};
  end

  I = bsxfun(@plus,h*(0:w-1),(1:h)');

  L = sparse( ...
    [vec(I(2:end,:)) vec(I(:,2:end))], ...
    [vec(I(1:end-1,:)) vec(I(:,1:end-1))], ...
    1,w*h,w*h);
  L = L + L';
  L = L-diag(sum(L,2));
end
