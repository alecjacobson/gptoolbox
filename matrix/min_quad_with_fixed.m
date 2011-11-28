function [Z,F] = min_quad_with_fixed(varargin)
  % MIN_QUAD_WITH_FIXED Minimize quadratic energy Z'*A*Z + Z'*B + C with
  % constraints that Z(known) = Y, optionally also subject to the constraints
  % Aeq*Z = Beq
  %
  % [Z,F] = min_quad_with_fixed(A,B,known,Y)
  % [Z,F] = min_quad_with_fixed(A,B,known,Y,Aeq,Beq)
  % [Z,F] = min_quad_with_fixed(A,B,known,Y,Aeq,Beq,F)
  %
  % Inputs:
  %   A  n by n matrix of quadratic coefficients
  %   B  n by 1 column of linear coefficients
  %   known  #known list of indices to known rows in Z
  %   Y  #known by cols list of fixed values corresponding to known rows in Z
  %   Optional:
  %     Aeq  m by n list of linear equality constraint coefficients
  %     Beq  m by 1 list of linear equality constraint constant values
  % Outputs:
  %   Z  n by cols solution
  %   Optional:
  %     F  struct containing all information necessary to solve a prefactored
  %     system touching only B, Y, and optionally Beq

  % Implementation details:
  % minimize x'Ax + x'B
  % subject to Aeq x = Beq
  %
  % This is the same as:
  % find the saddle point of x'Ax + x'B + lambda' * (Aeq * x  - Beq)
  % where lambda is a vector of lagrange multipliers, one for each of the
  % constraints (rows in Aeq)
  %
  % Then we rewrite this, combining x and lambda:
  % [x; lambda]' * [A Aeq';Aeq Z] * [x; lambda] + [x; lambda]' * [B; -2*Beq]
  %
  % Notice the -2 because lamba' * Aeq * x shows up twice in the quadratic part.
  % Then I can differentiate with respect to [x; lambda] and we get:
  % 2*[A Aeq';Aeq Z] * [x; lambda] + [B; -2*Beq]
  
  % Setting that to zero and moving the knowns to the right hand side we get:
  % [A Aeq';Aeq Z] * [x; lambda] = -0.5 * [B; -2*Beq]
  
  % Process input
  A = varargin{1};
  B = varargin{2};
  known = varargin{3};
  Y = varargin{4};
  Aeq = [];
  Beq = [];
  if nargin >= 6
    Aeq = varargin{5};
    Beq = varargin{6};
  end
  F = [];
  if nargin >= 7
    F = varargin{7};
  end

  if isempty(F)
    F = precompute(A,known,Aeq);
  end
  Z = solve(F,B,Y,Beq);

  function F = precompute(A,known,Aeq)
    % PRECOMPUTE perform any necessary precomputation of system including
    % factorization and preparation of right-hand side
    % 
    % F = precompute(A,known,Aeq)
    %
    % Inputs:
    %   A  n by n matrix of quadratic coefficients
    %   known  #known list of indices to known rows in Z
    %   Optional:
    %     Aeq  m by n list of linear equality constraint coefficients
    % Outputs:
    %   F  struct containing all information necessary to solve a prefactored
    %   system touching only B, Y, and optionally Beq

    % number of rows
    n = size(A,1);
    % cache problem size
    F.n = n;

    if isempty(Aeq)
      Aeq = zeros(0,n);
    end

    assert(size(A,1) == n);
    assert(size(A,2) == n);
    assert(prod(size(known)) == numel(known));
    assert(  isempty(known) || min(known) >= 1);
    assert(  isempty(known) || max(known) <= n);
    assert(n == size(Aeq,2));

    % cache known
    F.known = known;
    indices = 1:n;
    % get list of unknown variables including lagrange multipliers
    F.unknown = indices(~ismember(indices,known));

    Auu = A(F.unknown,F.unknown);

    % determine if A(unknown,unknown) is symmetric and/or postive definite
    %F.Auu_sym = ~any(any(Auu - Auu'));
    sym_measure = max(max(abs(Auu - Auu')))/max(max(abs(Auu)));
    %sym_measure = normest(Auu-Auu')./normest(Auu);
    if sym_measure > eps
      % not very symmetric
      F.Auu_sym = false;
    elseif sym_measure > 0
      % nearly symmetric but not perfectly
      F.Auu_sym = true;
    else
      assert(sym_measure == 0);
      % Perfectly symmetric
      F.Auu_sym = true;
    end
    % Determine if positive definite (also compute cholesky decomposition if it
    % is as a side effect)
    F.Auu_pd = false;
    if F.Auu_sym
      [F.L,p] = chol(Auu,'lower');
      F.Auu_pd = p==0;
    end
    F.Auu_sym
    F.Auu_pd

    % keep track of whether original A was sparse
    A_sparse = issparse(A);

    % check if there are blank constraints
    F.blank_eq = all(Aeq(:,F.unknown)==0,2);
    if any(F.blank_eq)
      warning('min_quad_with_fixed:blank_eq', [ ...
        'Removing blank constraints. ' ...
        'You ought to verify that known values satisfy contsraints']);
      Aeq = Aeq(~F.blank_eq,:);
    end
    % number of linear equality constraints
    neq = size(Aeq,1);
    assert(neq <= n);
    % get list of lagrange multiplier indices
    F.lagrange = n+(1:neq);

    Z = zeros(neq,neq);
    % append lagrange multiplier quadratic terms
    A = [A Aeq';Aeq Z];
    assert(A_sparse == issparse(A));
    % precompute RHS builders
    F.preY = A([F.unknown F.lagrange],known) + ...
      A(known,[F.unknown F.lagrange])';

    % create factorization
    if F.Auu_sym && F.Auu_pd
      if neq == 0
        % we already have F.L
        F.U = F.L';
      else
        % check that Auu is not ill-conditioned, if it is then don't us
        % lu_lagrange trick, rather use straight LU decomposition
        if log10(condest(Auu)) < 16
          [F.L,F.U,p] = lu_lagrange(Auu,Aeq(:,F.unknown)',F.L);
        else
          NA = A([F.unknown F.lagrange],[F.unknown F.lagrange]);
          [F.L,F.U] = lu(NA);
        end
      end
    else
      NA = A([F.unknown F.lagrange],[F.unknown F.lagrange]);
      % LU factorization of NA
      [F.L,F.U] = lu(NA);
    end
  end

  function Z = solve(F,B,Y,Beq)
    % SOLVE  perform solve using precomputation and parameters for building
    % right-hand side that are allowed to change without changing precomputation
    %
    % Z = solve(F,B,Y,Beq)
    %
    % Inputs:
    %   F  struct containing all information necessary to solve a prefactored
    %     system touching only B, Y, and optionally Beq
    %   B  n by 1 column of linear coefficients
    %   Y  #known by cols list of fixed values corresponding to known rows in Z
    %   Optional:
    %     Beq  m by 1 list of linear equality constraint constant values
    % Outputs:
    %   Z  n by cols solution
    %

    % number of known rows
    kr = numel(F.known);
    if kr == 0
      assert(isempty(Y));
      % force Y to have 1 column even if empty
      Y = zeros(0,1);
    end
    assert(kr == size(Y,1));
    cols = size(Y,2);

    if any(F.blank_eq)
      Beq = Beq(~F.blank_eq,:);
    end

    % number of lagrange multipliers aka linear equality constraints
    neq = numel(F.lagrange);

    if neq == 0
      assert(isempty(Beq));
      Beq = zeros(0,1);
    end

    % append lagrange multiplier rhs's
    B = [B; -2*Beq];

    NB = F.preY * Y + repmat(B([F.unknown F.lagrange]),1,cols);

    Z = zeros(F.n,cols);
    Z(F.known,:) = Y;

    Z([F.unknown F.lagrange],:) = -0.5 * (F.U \ (F.L \ NB)) ;

    if neq ~= 0
      % throw away lagrange multipliers
      Z = Z(1:(end-neq),:);
    end
    
  end
end
