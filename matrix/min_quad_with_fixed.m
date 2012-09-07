function [Z,F,Lambda,Lambda_known] = min_quad_with_fixed(varargin)
  % MIN_QUAD_WITH_FIXED Minimize quadratic energy Z'*A*Z + Z'*B + C with
  % constraints that Z(known) = Y, optionally also subject to the constraints
  % Aeq*Z = Beq
  % http://www.alecjacobson.com/weblog/?p=1913
  % http://www.alecjacobson.com/weblog/?p=2491
  %
  % [Z,F] = min_quad_with_fixed(A,B,known,Y)
  % [Z,F] = min_quad_with_fixed(A,B,known,Y,Aeq,Beq)
  % [Z,F] = min_quad_with_fixed(A,B,known,Y,Aeq,Beq,F)
  %
  % Inputs:
  %   A  n by n matrix of quadratic coefficients
  %   B  n by 1 column of linear coefficients, if empty then assumed B = 0
  %   known  #known list of indices to known rows in Z
  %   Y  #known by cols list of fixed values corresponding to known rows in Z
  %   Optional:
  %     Aeq  m by n list of linear equality constraint coefficients
  %     Beq  m by 1 list of linear equality constraint constant values, if
  %       empty then assumed Beq = 0
  % Outputs:
  %   Z  n by cols solution
  %   Optional:
  %     F  struct containing all information necessary to solve a prefactored
  %     system touching only B, Y, and optionally Beq
  %     Lambda  m list of values of lagrange multipliers corresponding to each
  %       row in Aeq
  %     Lambda_known  m list of values of lagrange multipliers corresponding to each
  %       known value
  %
  % Troubleshooting:
  %   'Warning: Matrix is singular to working precision' A number of things can
  %   cause this to happen:
  %     (1) Your system matrix A after removing knowns is not full rank, check
  %     condest of A and rank of A
  %     (2) Some constraints in Aeq are linearly dependent (after removing
  %     known values), check SVD of Aeq
  %   "My should be symmetric, positive definite but I see that lu is being
  %   used". Check (1) and (2) above, and:
  %     (3) be sure you giving A and not accidentally -A, which would be
  %     negative definite.
  %

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
  % treat empty B as column of zeros to match A
  if isempty(B)
    B = zeros(size(A,1),1);
  end
  known = varargin{3};
  Y = varargin{4};
  Aeq = [];
  Beq = [];
  if nargin >= 6
    Aeq = varargin{5};
    Beq = varargin{6};
  end
  % treat empty Beq as column of zeros to match Aeq
  if isempty(Beq)
    Beq = zeros(size(Aeq,1),1);
  end
  F = [];
  if nargin >= 7
    F = varargin{7};
  end

  if isempty(F)
    F = precompute(A,known,Aeq);
  end
  [Z,Lambda,Lambda_known] = solve(F,B,Y,Beq);

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

    assert(size(A,1) == n, ...
      'Rows of system matrix (%d) != problem size (%d)',size(A,1),n);
    assert(size(A,2) == n, ...
      'Columns of system matrix (%d) != problem size (%d)',size(A,2),n);
    assert(isempty(known) || min(size(known))==1, ...
      'known indices (size: %d %d) not a 1D list',size(known));
    assert(isempty(known) || min(known) >= 1, ...
      'known indices (%d) < 1',min(known));
    assert(isempty(known) || max(known) <= n, ...
      'known indices (%d) > problem size (%d)',max(known),n);
    assert(n == size(Aeq,2), ...
      'Columns of linear constraints (%d) != problem size (%d)',size(Aeq,2),n);

    % cache known
    F.known = known;
    % get list of unknown variables including lagrange multipliers
    %indices = 1:n;
    %F.unknown = indices(~ismember(indices,known));
    F.unknown = find(~sparse(1,known,true,1,n));

    Auu = A(F.unknown,F.unknown);
    % note that columns are in *original* order
    F.Ak = A(F.known,:);

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
    F.Auu_sym = false;

    % check if there are blank constraints
    F.blank_eq = ~any(Aeq(:,F.unknown),2);
    if any(F.blank_eq)
      warning('min_quad_with_fixed:blank_eq', [ ...
        'Removing blank constraints. ' ...
        'You ought to verify that known values satisfy contsraints']);
      Aeq = Aeq(~F.blank_eq,:);
    end
    % number of linear equality constraints
    neq = size(Aeq,1);
    %assert(neq <= n,'Number of constraints (%d) > problem size (%d)',neq,n);

    % Determine if positive definite (also compute cholesky decomposition if it
    % is as a side effect)
    F.Auu_pd = false;
    if F.Auu_sym && neq == 0
      % F.S'*Auu*F.S = F.L*F.L'
      [F.L,p,F.S] = chol(Auu,'lower');
      F.Auu_pd = p==0;
    end

    % keep track of whether original A was sparse
    A_sparse = issparse(A);

    % get list of lagrange multiplier indices
    F.lagrange = n+(1:neq);

    if neq > 0
      if issparse(A) && ~issparse(Aeq)
        warning('min_quad_with_fixed:sparse_system_dense_constraints', ...
        'System is sparse but constraints are not, solve will be dense');
      end
      if issparse(Aeq) && ~issparse(A)
        warning('min_quad_with_fixed:dense_system_sparse_constraints', ...
        'Constraints are sparse but system is not, solve will be dense');
      end
      Z = sparse(neq,neq);
      % append lagrange multiplier quadratic terms
      A = [A Aeq';Aeq Z];
      %assert(~issparse(Aeq) || A_sparse == issparse(A));
    end
    % precompute RHS builders
    F.preY = A([F.unknown F.lagrange],known) + ...
      A(known,[F.unknown F.lagrange])';

    % LDL has a different solve prototype
    F.ldl = false;
    % create factorization
    if F.Auu_sym
      if neq == 0 && F.Auu_pd
        % we already have F.L
        F.U = F.L';
        F.P = F.S';
        F.Q = F.S;
      else
      %    % check that Auu is not ill-conditioned, if it is then don't us
      %    % lu_lagrange trick, rather use straight LU decomposition
      %    if log10(condest(Auu)) < 16
      %      %% p will only be 0 if this works
      %      %[F.L,F.U,p] = lu_lagrange(Auu,Aeq(:,F.unknown)',F.L);
      %      %F.Q = 1;
      %      %F.P = F.Q';
      %      % 4 times faster (but only for small number of constraints in Aeq)
      %      [F.L,F.U,p,F.Q] = lu_lagrange(Auu,Aeq(:,F.unknown)',F.L,F.S);
      %      F.P = F.Q';
      %    else
      %      p = -1;
      %    end
      %    if p ~= 0
      %      NA = A([F.unknown F.lagrange],[F.unknown F.lagrange]);
      %      [F.L,F.U,F.P,F.Q] = lu(NA);
      %    end
        % LU_LAGRANGE is faster only for #constraints << #unknowns
        % LDL is faster than LU for moderate #constraints < #unknowns
        NA = A([F.unknown F.lagrange],[F.unknown F.lagrange]);
        [F.L,F.D,F.P,F.S] = ldl(NA);
        F.ldl = true;
      end
    else
      NA = A([F.unknown F.lagrange],[F.unknown F.lagrange]);
      % LU factorization of NA
      %[F.L,F.U] = lu(NA);
      %F.P = 1;
      %F.Q = 1;
      % 10 times faster
      if issparse(NA)
        [F.L,F.U,F.P,F.Q] = lu(NA);
      else
        [F.L,F.U] = lu(NA);
        F.P = 1;
        F.Q = 1;
      end
    end
  end

  function [Z,Lambda,Lambda_known] = solve(F,B,Y,Beq)
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
    %   Lambda  m by cols list of lagrange multiplier *values*
    %   Lambda_known  #known by cols list of lagrange multiplier *values* for
    %     known variables
    %

    % number of known rows
    kr = numel(F.known);
    if kr == 0
      assert(isempty(Y),'Known values should not be empty');
      % force Y to have 1 column even if empty
      if size(Y,2) == 0
        Y = zeros(0,1);
      end
    end
    assert(kr == size(Y,1), ...
      'Number of knowns (%d) != rows in known values (%d)',kr, size(Y,1));
    cols = size(Y,2);

    if any(F.blank_eq)
      Beq = Beq(~F.blank_eq,:);
    end

    % number of lagrange multipliers aka linear equality constraints
    neq = numel(F.lagrange);

    if neq == 0
      assert(isempty(Beq),'Constraint right-hand sides should not be empty');
      Beq = zeros(0,1);
    end

    %if cols ~= size(B,2)
    %  B = repmat(B,1,cols);
    %end
    %if cols ~= size(Beq,2)
    %  Beq = repmat(Beq,1,cols);
    %end
    %% append lagrange multiplier rhs's
    %B = [B; -2*Beq];
    %assert(size(Y,2) == size(B,2));
    %NB = F.preY * Y + B([F.unknown F.lagrange],:);
    %NBold = bsxfun(@plus, F.preY * Y, [B(F.unknown,:); -2*Beq(F.lagrange-F.n,:)]);
    NB = ...
      bsxfun(@plus, ...
        bsxfun(@plus,  ...
          F.preY * Y,  ...
          [B(F.unknown,:); zeros(numel(F.lagrange),size(B,2))]), ...
        [zeros(numel(F.unknown),size(Beq,2)); -2*Beq(F.lagrange-F.n,:)]);

    Z = zeros(F.n,cols);
    Z(F.known,:) = Y;

    if F.ldl
      Z([F.unknown F.lagrange],:) = ...
        -0.5 * F.S * (F.P * (F.L'\(F.D\(F.L\(F.P' * (F.S * NB))))));
    else
      Z([F.unknown F.lagrange],:) = -0.5 * F.Q * (F.U \ (F.L \ ( F.P * NB)));
    end

    Lambda = [];
    if neq ~= 0
      % save lagrange multipliers
      Lambda = Z(F.lagrange,:);
      % throw away lagrange multipliers
      Z = Z(1:(end-neq),:);
    end

    Lambda_known = -bsxfun(@plus,F.Ak * Z,0.5*B(F.known,:));
    
  end
end
