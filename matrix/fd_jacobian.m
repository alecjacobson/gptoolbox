function dfdx = fd_jacobian(obj,x,epsilon)
  % dfdx = fd_jacobian(obj,x,epsilon)
  %
  % Inputs:
  %  obj  function handle callable with f = obj(x)
  %  x  #x query point 
  % Outputs:
  %  dfdx  #f by #x finite difference approximation of the jacobian of obj at x
  % 
  % Example:
  %   x = [1;1;1];
  %   A = magic(size(x,1));
  %   A = A+A';
  %   obj = @(x) 0.5*x.'*A*x;
  %   % Hessian
  %   H = fd_jacobian(@(x) fd_jacobian(obj,x),x);
  %   max(abs(H(:)-A(:)))
  %
  %
  % Obsolete, use fd
  %
  % See also: fd
  x = reshape(x,[],1);
  if nargin<3
    epsilon = 1e-5;
  end
  out = obj(x);
  nout = numel(out);
  if issparse(out)
    dfdx = sparse(nout,numel(x));
  else
    dfdx = zeros(nout,numel(x));
  end
  for i = 1:numel(x)
    dx = zeros(numel(x),1);
    dx(i) = 1*epsilon;
    fp = obj(x + dx);
    fm = obj(x - dx);
    dfdx(:,i) = (fp-fm)/(2*epsilon);
  end
end

