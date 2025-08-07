function dfdx = fd(f,x,epsilon)
  % FD Finite difference approximation of the derivative of a function f
  %
  % dfdx = fd(f,x,epsilon)
  %
  % Inputs:
  %   f  function handle taking x as input and outputting a vector of size #f
  %   x  #x₁ by #x₂ by ... by #xₙ vector of inputs, where #x₁ is the first dimension
  %   epsilon  (optional) step size for finite difference, default 1e-6
  % Outputs:
  %   dfdx  #f by #x₁ by #x₂ by ... by #xₙ matrix of derivatives
  %
  % See also: fd_jacobian
  if nargin < 3
    epsilon = 1e-6;
  end
  % dfdx is [size(f) numel(x)]
  dfdx = [];
  for i = 1:numel(x)
    x_i = x(i);
    x(i) = x_i + epsilon;
    fxi_p = f(x);
    x(i) = x_i - epsilon;
    fxi_m = f(x);
    % restore
    x(i) = x_i;

    if i == 1
      size_f = size(fxi_p);
      xdim = find(size(fxi_p)~=1,1,'last');
      if isempty(xdim)
        xdim = 2;
      else
        xdim = xdim + 1;
      end
      size_f(xdim) = numel(x);
      dfdx = zeros(size_f);
    end


    dfdx_i = (fxi_p - fxi_m) / (2*epsilon);
    switch xdim
    case {1,2}
      dfdx(:,i) = dfdx_i;
    case 3
      dfdx(:,:,i) = dfdx_i;
    case 4
      dfdx(:,:,:,i) = dfdx_i;
    otherwise
      error('not implemented');
    end

  end


end
