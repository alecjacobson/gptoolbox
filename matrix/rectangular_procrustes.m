function X = rectangular_procrustes(A,B)
  % Solve a "procrustes" problem where the unkown matrix is not a sqaure matrix,
  % but rather has many more rows than columns. Specifically, we find a solution
  % to:
  %
  % min ‖A*X - B‖²  subject to X'*X = I
  %  X 
  %
  % Inputs:
  %   A  #A by #X (sparse) constraint matrix 
  %   B  #A by #B  matrix of right-hand sides
  % Outputs:
  %   X  #X by #B solution matrix
  %
  % Note: unlike the builtin `procrustes` function this does not allow for
  % translation or scaling, but does allow for "reflection". That is, we're
  % really just ensuring that X'*X = I, period.
  %

  % Inverse square root of a diagonal matrix
  isqrt = @(D) diag(sqrt(diag(D)).^-1);
  Q = A'*B;
  % Proof that least squares problem is minimized (though a direct reading would
  % lead to an expensive svd of a #A by #B matrix)
  % http://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec21.pdf
  [U,S,V] = svd( Q'*Q );
  X = A'*(B*(U*isqrt(S)*V'));
  % Proof that constraint is satisfied.
  % X = A'*B*U*√S⁻¹*V'
  % X'*X = (A'*B*U*√S⁻¹*V')' * (A'*B*U*√S⁻¹*V')
  % X'*X = V*√S⁻¹*U'* (A'*B)'*(A'*B) *U*√S⁻¹*V'
  % X'*X = V*√S⁻¹*U'*    U*S*V'      *U*√S⁻¹*V'
  %   Q'*Q is symmetric therefore U = V'
  %   TODO: This implies we should have just called eig.
  % X'*X = V*√S⁻¹*V*V'*S*V'*V'*√S⁻¹*V'
  % X'*X = V*√S⁻¹*S*√S⁻¹*V'
  % X'*X = V*V'
  % X'*X = I

end
