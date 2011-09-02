function [Z,F,sym] = min_quad_with_fixed(A,B,known,Y,Aeq,Beq,F,sym)
  % MIN_QUAD_WITH_FIXED Minimize quadratic energy Z'*A*Z + Z'*B + C with
  % constraints that Z(known) = Y, optionally also subject to the constraints
  % Aeq*Z = Beq
  %
  % [Z,F,sym] = min_quad_with_fixed(A,B,known,Y)
  % [Z,F,sym] = min_quad_with_fixed(A,B,known,Y,Aeq,Beq)
  % [Z,F,sym] = min_quad_with_fixed(A,B,known,Y,Aeq,Beq,F,sym)
  %
  % Inputs:
  %   A  n by n matrix of quadratic coefficients
  %   B  n by 1 column of linear coefficients
  %   known list of indices to known rows in Z
  %   Y  list of fixed values corresponding to known rows in Z
  %   Optional:
  %     Aeq  m by n list of linear equality constraint coefficients
  %     Beq  m by 1 list of linear equality constraint constant values
  %     F  struct containing factorization of system (from previous run of this
  %       function with identical A,Aeq matrices and known indices
  %     sym  flag specifying whether final system is symmetric
  % Outputs:
  %   Z  n by cols solution
  %   F  struct containing factorization of given problem's system matrix
  %   sym  flag specifying whether final system is symmetric
  %

  % number of rows
  n = size(A,1);
  % number of columns
  if(isempty(Y))
    cols = 1;
  else
    cols = size(Y,2);
  end

  assert(size(A,1) == n);
  assert(size(A,2) == n);
  assert(size(B,1) == n);
  assert(size(B,2) == 1);
  assert(prod(size(known)) == numel(known));
  % number of known rows
  kr = numel(known);
  assert(isempty(known) || min(known) >= 1);
  assert(isempty(known) || max(known) <= n);
  assert(kr == size(Y,1));

  % default is to have 0 linear equality constraints
  neq = 0;
  if(exist('Aeq','var') && ~isempty(Aeq))
    % number of linear equality constraints
    neq = size(Aeq,1);
    assert(size(Aeq,2) == n);
    assert(size(Beq,1) == neq);
    assert(size(Beq,2) == 1);
    assert(neq <= n);
    assert(issparse(A));
    assert(issparse(Aeq));
    % append linear constraints as lagrange multipliers
    [AI,AJ,AV] = find(A);
    [AeqI,AeqJ,AeqV] = find(Aeq);
    % conver to column vectors
    AeqI = AeqI(:);
    AeqJ = AeqJ(:);
    AeqV = AeqV(:);
    A = sparse([AI;AeqI+n;AeqJ], [AJ;AeqJ;AeqI+n],[AV;AeqV;AeqV], n+neq,n+neq);
    B = [B; -2*Beq];
    % set number of variables to include lagrange multipliers
    n = n+neq;
  end


  % determine symmetry of final system
  if ~exist('sym','var') || isempty(sym)
    sym = ~any(any(A - A'));
  end

  indices = 1:n;
  % get list of unknown variables
  unknown = indices(~ismember(indices,known));

  % new quadratic and linear terms of energy without known values (can be
  % improved for symmetric matrices since A(unknown,known) == A(known,unknown)'
  NA = A(unknown,unknown);

  % can it be that A is symmetric but NA is not... certainly
  assert( sym ==  ~any(any(NA - NA')));

  if(isempty(known))
    NB = repmat(B(unknown),1,cols);
  else
    if sym
      NB = 2*A(unknown,known) * Y + repmat(B(unknown),1,cols);
    else
      NB = A(unknown,known) * Y + A(known,unknown)' * Y + ...
        repmat(B(unknown),1,cols);
    end
  end

  Z = zeros(n,cols);
  Z(known,:) = Y;

  % factor NA
  if ~exist('F','var') || isempty(F)
    if sym
      [F.L, g, F.PT] = chol(NA, 'lower');
      % be sure chol factorization worked
      if(g ~= 0)
        %% superfluous warning
        %warning(['Resulting system was not positive definite, ' ...
        %  'trying negative system...']);
        %% if it didn't then NA is not positive definite, maybe it's
        %% "negative-definite" so try letting NA=-NA (thus NB = -NB)
        %NA = -NA;
        %NB = -NB;
        %[F.L, g, F.PT] = chol(NA, 'lower');
        %assert(g == 0);
        error(['Resulting system is not positive-definite ' ... 
          '(perhaps you gave -A instead of A)']);
      end
    else
      [F.L,F.U,F.P,F.Q,F.R] = lu(NA);
    end
  end

  % solve 
  if sym
    Z(unknown,:) = -0.5 * F.PT * (F.L' \ (F.L \ (F.PT' * NB))) ;
  else
    Z(unknown,:) = -0.5*(F.Q*(F.U\(F.L\(F.P*(F.R\NB)))));
  end
  % Z(unknown,:) = -0.5*(NA\NB);


  if neq ~= 0
    % throw away lagrange multipliers
    Z = Z(1:(n-neq),:);
  end
  
end
