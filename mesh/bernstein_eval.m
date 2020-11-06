function B = bernstein_eval(i,n,u)
  % BERNSTERIN_EVAL evaluate a Bernstein polynomial. B_i,n(u)
  %
  % B = bernstein_eval(i,n,u)
  % 
  % Inputs:
  %  i  #i scalar index of polynomial
  %  n  scalar order of polynomial
  %  u  #u evaluation points
  % Outputs:
  %  B  #u by #i values
  %   
  u = reshape(u,numel(u),1);
  i = reshape(i,1,numel(i));
  nck = factorial(n)./(factorial(i).*factorial(n-(i)));
  B = nck .*  u.^i .* (1-u).^(n-i);
end
