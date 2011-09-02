function A,rhs = curve_editing_1d_system( undeformed, h)
  % Constructs the system matrix and right-hand side for minimizing the
  % enrgy function:
  % F (y) = (1/2)\int{y"(t)}{dt}
  % which is to say it solves
  % y???? = 0
  %
  % The approximation for the fourth derivative used is:
  % (fj-2 - 4*fj-1 + 6*fj - 4*fj+1 + fj+2)/(h^4)
  %
  % NOTE: the returned system is missing boundary and fixed point
  % conditions
  %
  % INPUT
  % undeformed  original input positions
  % h           distance between elements
  %
  % OUTPUT
  % A           system matrix
  % rhs         right-hand side
  
  N = size(undeformed,1);
  %A = hessian_direct(N);
  A = zeros(N,N);
  enum = 1:N;
  A(enum,enum) = 
  coeff = 1./(h*h*h*h);
  
  rhs = zeros([N, 1]);
  
  
  
end