function [t,x,fx] = backtracking_line_search(f,x0,dfx0,dx,alpha,beta,max_iter)
  % [t,x,fx] = backtracking_line_search(f,x0,dfx0,dx,alpha,beta,max_iter)
  %
  % Inputs:
  %   f  function such that f(x) returns scalar objective to minimize
  %   x0  #x initial solution
  %   dfx0  gradient of f at x0
  %   dx  #x proposed step in x
  %   alpha  the fraction of the decrease in f predicted by linear extrapolation
  %      that we will accept. Must be ∈(0,½), usually between 0.01 and 0.3.
  %   beta  decrease factor per iteration ∈(0,1), usually between 0.1 and 0.8.
  %   max_iter  maximum number of iterations {30}
  % Outputs:
  %   t  selected step size
  %   x  #x vector of solution after step: x = x0+t*dx
  %   fx  scalar objective after step fx = f(x)
  %

  % [Boyd 2006] "Convex Optimization" Chapter 11 "Interior-point methods"
  % Algorithm 9.2
  assert(alpha>0 && alpha<0.5);
  assert(beta>0 && beta<1);
  t = 1;
  fx0 = f(x0);
  if nargin<7
    max_iter = 30;
  end
  for iter = 1:max_iter
    x = x0+t*dx;
    fx = f(x);
    if fx<=fx0+alpha*t*dfx0'*dx
      return;
    end
    t = beta*t;
  end
  t = 0;
  x = x0;
  fx = fx0;
end
