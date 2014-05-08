function [Q,J,F,T,I,V,E] = spinv(A,varargin)
  % SPINV Compute a sparse factorization of the Moore-Penrose pseudoinverse A⁺
  % of a sparse matrix A.  This is useful for solving
  % under-and-over-constrained systems Ax = b: where A is m by n with
  % rank r < min(m,n). Then to solve the system:
  %
  %     x = E*(V*(I'*(T'\(F'*(J'*(Q'*b))))))
  %    
  % Note: Store these sparse factors and sovle using the routine above **with
  % parentheses**. Do not precompute A⁺ as it will be dense.
  %
  % Inputs:
  %   A  m by n sparse
  %   Optional:
  %     'Epsilon' followed by epsilon used to determine rank from qr {1e-15}
  % Outputs:
  %   Q  m by m sparse orthogonal matrix
  %   J  m by r sparse rectangular identity matrix
  %   F  r by r sparse permutation matrix
  %   T  r by r sparse triangular matrix (diagonal?)
  %   I  r by n sparse rectangular identity matrix
  %   V  n by n sparse orthogonal matrix
  %   E  n by n sparse permutation matrix
  %   
  % Example:
  %   % over and under determined system of 3 variables
  %   A = [1 0 0;1 0 0;0 1 1;0 1 1];
  %   b = [4;5;6;7];
  %   % Compute facotrization of sparse pseudo inverse
  %   [Q,J,F,T,I,V,E] = spinv(A);
  %   % Compute minimum norm best fit
  %   xs = E*(V*(I'*(T'\(F'*(J'*(Q'*b))))));
  %   % Compare to Tikhonov regularization
  %   xt = (A'*A+1e-10*speye(3,3))\(A'*b);
  %   

  % default values
  epsilon = 1e-15;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'Epsilon'}, {'epsilon'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace 
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end
  % http://www.irisa.fr/sage/wg-statlin/WORKSHOPS/LEMASSOL05/SLIDES/QR/Guyomarch.pdf

  % dimensions of A
  m = size(A,1);
  n = size(A,2);
  % QR factorization of AE
  [Q,R,E] = qr(A);
  % We have that:
  %     AE = QR
  % rank
  r = find(any(abs(R)>epsilon,2),1,'last');
  % Only keep upper, non-zero block
  R1R2 = R(1:r,:);
  J = speye(m,r);
  % We have that:
  %    AE = Q*J*R1R2
  % QR factorization of R1R2'
  [V,Tfull,F] = qr(R1R2');
  % This means:
  %    R1R2'F = V*Tfull
  % only keep upper, non-zero block
  T = Tfull(1:r,:);
  I = speye(r,n);
  % and now:
  %    R1R2'F = V*I'*T
  %    R1R2' = V*I'*T*F'
  % Now we have that
  %     AE = Q*J*F*T'*I*V'
  % We want to solve:
  %     Ax = b
  %     AEE'x = b
  %     E'x = (Q*J*F*T'*I*V')⁺b
  %     x = E*(Q*J*F*T'*I*V')⁺b
  %     x = E*(V*(I'*(T'\(F'*(J'*(Q'*b))))))
  % 
end
