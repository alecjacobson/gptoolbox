function dfdx = cs_jacobian(obj,x,epsilon)
  % CS_JACOBIAN Compute the jacobian of a real-valued function which is
  % holomorphic near x using complex step finite-differences.
  % (easy to not be holomorphic accidentally:
  % https://www.alecjacobson.com/weblog/4921.html)
  %
  % dfdx = cs_jacobian(obj,x,epsilon)
  %
  % Inputs:
  %  obj  function handle callable with f = obj(x)
  %  x  #x query point 
  % Outputs:
  %  dfdx  #f by #x finite difference approximation of the jacobian of obj at x
  % 
  %   x = [1;1;1];
  %   A = magic(size(x,1));
  %   A = A+A.';
  %   % using `.'` instead of `'` is super important here!
  %   obj = @(x) 0.5*x.'*A*x;
  %   % Hessian
  %   H = fd_jacobian(@(x) cs_jacobian(obj,x),x);
  %   max(abs(H(:)-A(:)))
  %
  % Obsolete, use fd(...,use_complex_step=true)
  % 
  % See also: fd
  x = reshape(x,[],1);
  if nargin<3
    epsilon = 1e-100;
  end
  dfdx = zeros(numel(obj(x)),numel(x));
  for i = 1:numel(x)
    dx = zeros(numel(x),1);
    dx(i) = 1i*epsilon;
    fp = obj(x + dx);
    imag_fp = imag(fp);
    dfdx(:,i) = imag(fp)/epsilon;
  end
end

