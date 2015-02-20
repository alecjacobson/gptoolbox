function [Z,F,Lambda,Lambda_known] = min_quad_with_fixed(A,B,known,Y,Aeq,Beq,F)
  % MIN_QUAD_WITH_FIXED Minimize quadratic energy Z'*A*Z + Z'*B + C with
  % constraints that Z(known) = Y, optionally also subject to the constraints
  % Aeq*Z = Beq
  % http://www.alecjacobson.com/weblog/?p=1913
  % http://www.alecjacobson.com/weblog/?p=2491
  %
  % [Z,F] = min_quad_with_fixed(A,B,known,Y)
  % [Z,F,Lambda,Lambda_known] = min_quad_with_fixed(A,B,known,Y,Aeq,Beq,F)
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
  %     F see output
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
  %   "My hessian should be symmetric, positive definite but I see that lu is
  %     being used". Check (1) and (2) above, and:
  %     (3) be sure you giving A and not accidentally -A, which would be
  %     negative definite.
  %   "My constraints are not satisfied. That is, some abs(Aeq * Z - Beq) are
  %     not zero." In the output, check F.Aeq_li. If this is false then
  %     according to QR decomposition your constraints are linearly dependent.
  %     Check that your constraints are not conflicting.  Redundant or linearly
  %     dependent constraints **equations** (including rhs) should be OK, but
  %     linearly dependent rows in Aeq with mismatching rows in Beq mean
  %     there's a conflict.
  %     
  %
  %
  % Example:
  %   % one-linear to use pcg with same prototype:
  %   min_quad_with_fixed_pcg = @(A,B,known,Y,tol,iter,fun) ...
  %     full(sparse( ...
  %       [setdiff((1:size(A,1))',known(:));known(:)],1, ...
  %       [pcg( ...
  %         A(setdiff(1:end,known),setdiff(1:end,known)), ...
  %         -(A(setdiff(1:end,known),known) * Y) - ...
  %         0.5*B(setdiff(1:end,known)),tol,iter);bc]));

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
  
  if nargin < 4
    Y = [];
    known = [];
  end
  if nargin < 6
    Aeq = [];
    Beq = [];
  end
  % treat empty Beq as column of zeros to match Aeq
  if isempty(Beq)
    Beq = zeros(size(Aeq,1),1);
  end
  if nargin < 7
    F = [];
  end
  
  if isempty(F) || ~isfield(F,'precomputed') || F.precomputed == false
%     if ~isempty(F)
%       warning('Precomputing');
%     end
    F = precompute(A,known,Aeq,F);
  end
  [Z,Lambda,Lambda_known] = solve(F,B,Y,Beq);

  % !!SHOULD REMOVE F AS INPUT PARAM!!
  function F = precompute(A,known,Aeq,F)
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
    F.unknown = find(~sparse(1,known,true,1,n));

    Auu = A(F.unknown,F.unknown);
    % note that columns are in *original* order
    F.Ak = A(F.known,:);


    % determine if A(unknown,unknown) is symmetric and/or postive definite
    sym_measure = max(max(abs(Auu - Auu')))/max(max(abs(Auu)));
    %sym_measure = normest(Auu-Auu')./normest(Auu);
    if sym_measure > eps
      % not very symmetric
      F.Auu_sym = false;
    elseif sym_measure > 0
      % nearly symmetric but not perfectly
      F.Auu_sym = true;
    else
      % Either Auu is empty or sym_measure should be perfect
      assert(isempty(sym_measure) || sym_measure == 0,'not symmetric');
      % Perfectly symmetric
      F.Auu_sym = true;
    end

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
      if issparse(Auu)
        [F.L,p,F.S] = chol(Auu,'lower');
      else
        [F.L,p] = chol(Auu,'lower');
        F.S = eye(size(F.L));
      end
      F.Auu_pd = p==0;
    end

    % keep track of whether original A was sparse
    A_sparse = issparse(A);

    % Determine number of linearly independent constraints
    if neq > 0 && ~(isfield(F,'force_Aeq_li') && ~isempty(F.force_Aeq_li)&& F.force_Aeq_li)
      %tic;
      % Null space substitution with QR
      [AeqTQ,AeqTR,AeqTE] = qr(Aeq(:,F.unknown)');
      nc = find(any(AeqTR,2),1,'last');
      if isempty(nc)
        nc = 0;
      end
      %fprintf('QR: %g secs\n',toc);
      assert(nc<=neq);
      F.Aeq_li = nc == neq;
    else
      F.Aeq_li = true;
    end
    if neq > 0 && isfield(F,'force_Aeq_li') && ~isempty(F.force_Aeq_li)
      F.Aeq_li = F.force_Aeq_li;
    end

    % Use raw Lagrange Multiplier method only if rows of Aeq are Linearly
    % Independent
    if F.Aeq_li
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
          % LDL is faster than LU for moderate #constraints < #unknowns
          NA = A([F.unknown F.lagrange],[F.unknown F.lagrange]);
          assert(issparse(NA));
          [F.L,F.D,F.P,F.S] = ldl(NA);
          F.ldl = true;
        end
      else
        NA = A([F.unknown F.lagrange],[F.unknown F.lagrange]);
        % LU factorization of NA
        if issparse(NA)
          [F.L,F.U,F.P,F.Q] = lu(NA);
        else
          [F.L,F.U] = lu(NA);
          F.P = 1;
          F.Q = 1;
        end
      end
    else
      % We alread have CTQ,CTR,CTE
      %tic;
      % Aeq' * AeqTE = AeqTQ * AeqTR
      % AeqTE' * Aeq = AeqTR' * AeqTQ'
      % Aeq x = Beq
      % Aeq (Q2 lambda + lambda_0) = Beq
      % we know Aeq Q2 = 0 --> Aeq Q2 lambda = 0
      % Aeq lambda_0 = Beq
      % AeqTE' * Aeq lambda_0 = AeqTE' * Beq
      % AeqTR' * AeqTQ' lambda_0 = AeqTE' * Beq
      % AeqTQ' lambda_0 = AeqTR' \ (AeqTE' * Beq)
      % lambda_0 = AeqTQ * (AeqTR' \ (AeqTE' * Beq))
      % lambda_0 = Aeq \ Beq;
      % lambda_0 = AeqTQ * (AeqTR' \ (AeqTE' * Beq));
      AeqTQ1 = AeqTQ(:,1:nc);
      AeqTR1 = AeqTR(1:nc,:);
      %lambda_0 = [AeqTQ1 * (AeqTR1' \ (AeqTE' * Beq))];
      %fprintf('lambda_0: %g secs\n',toc);
      %tic;
      % Substitute x = Q2 lambda + lambda_0
      % min 0.5 x' A x - x' b
      %   results in A x = b
      % min 0.5 (Q2 lambda + lambda_0)' A (Q2 lambda + lambda_0) - (Q2 lambda + lambda_0)' b
      % min 0.5 lambda' Q2' A Q2 lambda + lambda Q2' A lambda_0 - lambda Q2' b 
      %  results in Q2' A Q2 lambda = - Q2' A lambda_0 + Q2' b
      AeqTQ2 = AeqTQ(:,(nc+1):end);
      QRAuu =  AeqTQ2' * Auu * AeqTQ2;
      %QRb = -AeqTQ2' * Auu * lambda_0 + AeqTQ2' * b;
      % precompute RHS builders
      F.preY = A(F.unknown,known) + A(known,F.unknown)';
      %fprintf('Proj: %g secs\n',toc);
      %tic;
      % QRA seems to be PSD
      [F.L,p,F.S] = chol(QRAuu,'lower');
      F.U = F.L';
      F.P = F.S';
      F.Q = F.S;
      %fprintf('Chol: %g secs\n',toc);
      % Perhaps if Auu is not PD then we need to use LDL...
      assert(p==0);
      % WHICH OF THESE ARE REALLY NECESSARY?
      F.Aeq = Aeq;
      F.AeqTQ2 = AeqTQ2;
      F.AeqTQ1 = AeqTQ1;
      F.AeqTR1 = AeqTR1;
      F.AeqTE = AeqTE;
      F.Auu = Auu;
    end
    F.precomputed = true;
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
    if isempty(Y)
      % use linear coefficients to determine cols
      if isempty(B)
        if isempty(Beq)
          cols = 1;
          Beq = zeros(0,cols);
        else
          cols = size(Beq,2);
        end
        B = zeros(F.n,cols);
      else
        cols = size(B,2);
      end
      Y = zeros(0,cols);
    else
      cols = size(Y,2);
      if isempty(B)
        B = zeros(F.n,cols);
      end
    end

    if any(F.blank_eq)
      Beq = Beq(~F.blank_eq,:);
    end


    % Build system's rhs
    if F.Aeq_li
      % number of lagrange multipliers aka linear equality constraints
      neq = numel(F.lagrange);
      if neq == 0
        assert(isempty(Beq),'Constraint right-hand sides should not be empty');
        Beq = zeros(0,1);
      end

      NB = ...
        bsxfun(@plus, ...
          bsxfun(@plus,  ...
            F.preY * Y,  ...
            [B(F.unknown,:); zeros(numel(F.lagrange),size(B,2))]), ...
          [zeros(numel(F.unknown),size(Beq,2)); -2*Beq(F.lagrange-F.n,:)]);
          

      % prepare solution
      Z = zeros(F.n+neq,cols);
      Z(F.known,:) = Y;

      if F.ldl
        Z([F.unknown F.lagrange],:) = ...
          -0.5 * F.S * (F.P * (F.L'\(F.D\(F.L\(F.P' * (F.S * NB))))));
      else
        Z([F.unknown F.lagrange],:) = -0.5 * F.Q * (F.U \ (F.L \ ( F.P * NB)));
      end

      % fix any removed constraints (set Lambda to 0)
      Lambda = zeros(numel(F.blank_eq),cols);
      if neq ~= 0
        % save lagrange multipliers
        Lambda(~F.blank_eq,:) = Z(F.lagrange,:);
        % throw away lagrange multipliers
        Z = Z(1:(end-neq),:);
      end
    else
      % Adjust Aeq rhs to include known parts
      Beq = -F.Aeq(:,known)*Y + Beq;
      % Where did this -0.5 come from? Probably the same place as above.
      NB = -0.5*(B(F.unknown,:) + F.preY * Y);
      eff_Beq = F.AeqTE' * Beq;
      % can't solve rectangular system: trim (expects that constraints are not
      % contradictory)
      AeqTR1T = F.AeqTR1';
      AeqTR1T = AeqTR1T(1:size(F.AeqTQ1,2),1:size(F.AeqTQ1,2));
      eff_Beq = eff_Beq(1:size(F.AeqTQ1,2));
      lambda_0 = F.AeqTQ1 * (AeqTR1T \ eff_Beq);
      QRB = -F.AeqTQ2' * (F.Auu * lambda_0) + F.AeqTQ2' * NB;
      lambda = F.Q * (F.U \ (F.L \ ( F.P * QRB)));
      % prepare solution
      Z = zeros(F.n,cols);
      Z(F.known,:) = Y;
      Z(F.unknown) = F.AeqTQ2 * lambda + lambda_0;
      Aequ = F.Aeq(:,F.unknown);
      % http://www.math.uh.edu/~rohop/fall_06/Chapter3.pdf
      %Lambda = (F.AeqTQ1' * Aequ') \ (F.AeqTQ1' * NB - F.AeqTQ1' * F.Auu * Z(F.unknown));
      % Can't solve rectangular system
      %Lambda = F.AeqTE * (F.AeqTR1 \ (F.AeqTQ1' * NB - F.AeqTQ1' * F.Auu * Z(F.unknown)));
      % TRIM: (other linearly dependent constraints get 0s?)
      Lambda = F.AeqTE * [ ...
        (F.AeqTR1(:,1:size(F.AeqTR1,1)) \ ...
          (F.AeqTQ1' * NB - F.AeqTQ1' * F.Auu * Z(F.unknown))); ...
        zeros(size(F.AeqTE,2)-size(F.AeqTR1,1),1)
        ];
    end

    Lambda_known = -bsxfun(@plus,F.Ak * Z,0.5*B(F.known,:));
  end

end
