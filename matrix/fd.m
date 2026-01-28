function dfdx = fd(f,x,epsilon,use_complex_step)
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
  % Example:
  %   % pointwise derivatives
  %   [X,Y] = meshgrid(linspace(-1,1),linspace(-1,1));
  %   Z = @(X,Y) 1-sqrt(X.^2 + Y.^2);
  %   dZdXY = fd(@(x) Z(X+x(1), Y+x(2)), [0 0]);
  %   surf(X,Y,Z(X,Y),'CData',atan2(dZdXY(:,:,1), dZdXY(:,:,2)));
  %   colormap(okloop(256));
  %
  % See also: fd_jacobian
  if nargin < 3
    epsilon = 1e-6;
  end
  if nargin < 4
    use_complex_step = false;
  end
  % dfdx is [size(f) numel(x)]
  dfdx = [];
  dx = zeros(size(x));
  for i = 1:numel(x)
    if use_complex_step
      dx(i) = 1i * epsilon;
      fxi_p = f(x + dx);
      imag_fp = imag(fxi_p);
      dfdx_i = imag_fp / epsilon;
      % restore
      dx(i) = 0;
    else
      x_i = x(i);
      x(i) = x_i + epsilon;
      fxi_p = f(x);
      x(i) = x_i - epsilon;
      fxi_m = f(x);
      dfdx_i = (fxi_p - fxi_m) / (2*epsilon);
      % restore
      x(i) = x_i;
    end

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
