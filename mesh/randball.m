function x = randball(n,dim)
  % RANDBALL Uniformly sample n points _inside_ a dim-dimensional ball.
  %
  % x = randball(n,dim)
  %
  % Inputs:
  %   n  number of samples
  %   dim  dimension of ball
  % Outputs:
  %   x  n by dim list of sample points
  % 
  x = normalizerow(normrnd(0,1,n,dim)).*(rand(n,1).^(1/dim));
end
