function y = callback_active_contour(x, options)

% callback_active_contour - callback for conjugate gradient
%   
%   y = callback_active_contour(x, options);
%
%   Copytight (c) 2007 Gabriel Peyre

% norm of gradient
dt = options.dt;
n = options.n;
d = options.d;
x = reshape(x,n,n);
y = divgrad( divgrad(x,options)./repmat(d,[1 1 2]),options );
if isempty(options.E)
    y = x(:) - dt * d(:) .* y(:);
else
    y = x(:) - dt * options.E(:) .* d(:) .* y(:);
end