function [varargout] = cubic_cubic_integrated_distance( ...
    x, ...
    y, ...
    g_C, ...
    h_C, ...
    P, ...
    g_D, ...
    h_D, ...
    Q)
  %
  % [E] = cubic_cubic_integrated_distance( x, y, g_C, h_C, P, g_D, h_D, Q)
  % [H,F,c] = cubic_cubic_integrated_distance( x, y, g_C, h_C, P, g_D, h_D)
  %
  % E = ½ ∫ₓʸ ‖ C( g_C u + h_C ) - D( g_D u + h_D ) ‖² du
  %
  % where C's control points are in rows of P and D's control points are in rows
  % of Q.
  %
  % E = 0.5*trace(Q.'*H*Q) + trace(Q.'*F) + c;
  % 
  % That is, 
  %   Q★ = argmin_Q E(Q) 
  %   Q★ = H⁻¹ F
  %   
  % Where
  %
  % C(T) = (1-T)³ Pᵢ⁰ + 3(1-T)²T Pᵢ¹ + 3(1-T)T² Pᵢ² + T³ Pᵢ³
  % D(T) = (1-T)³ Qᵢ⁰ + 3(1-T)²T Qᵢ¹ + 3(1-T)T² Qᵢ² + T³ Qᵢ³
  %
  % Inputs:
  %   x  starting value of integral 0≤x<1
  %   y  ending value of integral x<y≤1
  %   g_C  scalar multiplicative factor in parametric values of fixed curve
  %   h_C  scalar additive factor in parametric values of fixed curve
  %   P  4 by dim list of fixed curve control points
  %   g_D  scalar multiplicative factor in parametric values of fixed curve
  %   h_D  scalar additive factor in parametric values of fixed curve
  %   Q  4 by dim list of unknown curve control points (optional)
  % Outputs
  %   E  integrated energy
  %   H  4 by 4 quadratic form of energy as function of Q
  %   F  4 by dim linear term of energy as function of Q
  %   c  scalar constant of energy as function of Q
  %
  if isfloat(x) && isfloat(y)
    assert(x>=0);
    assert(x<1);
    assert(y>=x);
    assert(y<=1);
  end
  dim = size(P,2);

  % https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
  % n is number of points
  % 2n-1 order polynomial is integrated exactly
  % ↓
  % n = 4

  % https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules/quadrature_rules.html
  % n=4 Gauss-Legendre quadrature on [-1,1]
  w = [
    0.347854845137453857373063949222
    0.652145154862546142626936050778
    0.652145154862546142626936050778
    0.347854845137453857373063949222
    ];
  a = [
   -0.861136311594052575223946488893
   -0.339981043584856264802665759103
    0.339981043584856264802665759103
    0.861136311594052575223946488893
  ];
  % Move to [x,y] range.
  u = (0.5*a + 0.5)*(y-x)+x;
  % Adjust weights.
  w = ((y-x)/(1- -1))*w;

  % Evaluate C at quadrature points
  C = cubic_eval( P, g_C*u + h_C);


  % Map quadrature points to D's direct paramter space
  T = g_D * u + h_D;

  if nargout<=1
    % Compute energy directly. This is faster. Better be the same as below.
    % 
    % It's ever so slightly different (between 1e-16 and 1e-20)
    E = 0.5*sum(w.*(C - cubic_eval( Q, T)).^2,'all');
    varargout{1} = E;
    return;
  end
  
  % Stack bezier basis matrices for each evaluation point
  M = [(1-T).^3 3*T.*(1-T).^2 3*T.^2.*(1-T) T.^3];

  W = diag(w);

  H = M.'*W*M;
  assert(all(size(H) == [4 4]));
  F = -M.'*W*C;
  assert(size(F,1) == 4);
  assert(size(F,2) == dim);

  c = 0.5*trace(C.'*W*C);

  if isempty(Q)
    warning("Cant compute energy on empty Q");
    E = [];
  else
    % Actually compute energy
    E = 0.5*trace(Q.'*H*Q) + trace(Q.'*F) + c;
  end

  varargout{1} = H;
  varargout{2} = F;
  varargout{3} = c;
  varargout{4} = E;
end
