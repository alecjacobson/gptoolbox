function R = rand_rotation(n)
  % Uniformly sample SO(n) 
  %
  % R = rand_rotation(n)
  %
  % Inputs:
  %   n  dimension
  % Output:
  %   R  n by n rotation matrix
  %

  % This is 5× slower in 3D.
  %if n==3
  %  R = quat2mat(compact(randrot));
  %  return;
  %end

  % "How to generate a random unitary matrix" [Maris Ozols 2006]
  % http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
  [Q,R] = qr(randn(n));
  r = diag(R);
  L = diag(r./abs(r));
  R = Q*L;
  % I don't think this randomness is necessary
  %i = randperm(n,1);
  i = 1;
  R(i,:) = R(i,:)*det(R);
  assert(abs(det(R)-1)<1e-10)

end
