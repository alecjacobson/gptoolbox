function [x,err,it] = perform_conjugate_gradient(A,y,options)


% perform_conjugate_gradient - perform conjugate gradient
%
%     [x,err,k] = perform_conjugate_gradient(A,y,options);
%
%   Solves for A*x=y.
%   Works both for vector x,y and matrix x,y (parralel soving of many
%   linear equations).
%
%   Works only for semi-definite matrices.
%
%   A can be a matrix or a callback function y=A(x,options).
%   In this case, you have to set options.ncols as the number of columns of
%   A (if it is a callback).
%
%   err monitos the decreasing of |A*x-y|.
%   k is the total number of iterations.
%
%   You can set:
%       options.x is an initial guess
%       options.epsilon is maximum error
%       options.niter_max is maximum number of error
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
niter = getoptions(options, 'niter_max', 100);
epsilon = getoptions(options, 'epsilon', 1e-5);

use_callback = 0;
if not(isnumeric(A))
    use_callback = 1;
end
if isfield(options, 'x')
    x = options.x;
else
    if use_callback==0
        x = zeros(size(A,2),1);
    else
        if isfield(options, 'ncols')
            x = zeros(options.ncols, 1);
        else
            error('You have to specify options.ncols');
        end
    end
end


if norm(y, 'fro') == 0
    normb = epsilon;
else
    normb = epsilon * sum(y(:).^2);
end

if use_callback==0
    r  = y - A*x;
else
    r  = y - feval(A,x,options);    
end
r0 = sum(r.^2);
err = [sum(r0)];
for it=1:niter
    if (it==1)
        p = r;
    else
        % search direction
        beta = r0./rprev;
        p = r + repmat(beta, [size(x,1) 1]).*p;
    end
    % auxiliary vector
    if use_callback==0
        w  = A*p;
    else
        w  = feval(A,p,options);
    end
    d = sum(p .* w);
    I = find(abs(d)<eps); d(I) = 1;
    alpha = repmat( r0 ./ d, [size(x,1) 1] );              % found optimal alpha in line search
    x = x + alpha.*p;                       % new guess
    r = r - alpha.*w;                       % the residual is in fact r=b-A*x
    rprev = r0;                             % save norm of the old residual
    r0 = sum(r.^2);                         % compute norm of new residual
    err(end+1) = sum(r0);
    if err(end)<normb
        return;
    end
end



return;



err = [sqrt(r0)];
while ( sqrt(r0) > normb && k < kmax)
    if( k==1 ) 
        p = r;
    else
        beta = r0/rprev;
        % search direction
        p = r + beta*p;
    end
    if use_callback==1
        w = feval(A,p,options); 
    else
        w = A * p;                          %% auxiliary vector
    end
    alpha = r0 / (p' * w);              %% found optimal alpha in line search
    x = x + alpha*p;                    %% new guess
    r = r - alpha*w;                    %% the residual is in fact r=b-A*x
    rprev = r0;                         %% ... save norm of the old residual
    r0 = ( norm(r, 'fro') )^2;                 %% ... compute norm of new residual
    k = k + 1;    
    err(end+1) = sqrt(r0);
end