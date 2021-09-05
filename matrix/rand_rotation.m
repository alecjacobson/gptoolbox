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

  % "How to generate a random unitary matrix" [Maris Ozols 2006]
  % http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
  [Q,~] = qr(randn(n));
  s = (2*(rand>0.5)-1);
  Q(:,1) = Q(:,1)*s;
  %b = det(Q);
  b = s*(2*mod(n,2)-1);
  Q(:,2) = b*Q(:,2);
  assert(abs(det(Q)-1)<1e-10)
  % rename as R
  R = Q;
end
