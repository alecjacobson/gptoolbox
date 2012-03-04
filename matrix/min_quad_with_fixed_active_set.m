function [Z,F,Lambda,active_set_ieq,active_set_lx,active_set_ux] = min_quad_with_fixed_active_set(varargin)
  % MIN_QUAD_WITH_FIXED_ACTIVE_SET Minimize quadratic energy Z'*A*Z + Z'*B + C
  % with constraints that Z(known) = Y, optionally also subject to the
  % constraints Aeq*Z = Beq, and further optionally subject to the linear
  % inequality constraints that Aieq*Z <= Bieq and constant inequality
  % constraints lx <= x <= ux
  %
  % [Z] = min_quad_with_fixed_active_set(A,B,known,Y,Aeq,Beq,Aeiq,Beiq,lx,ux)
  %
  % Inputs:
  %   A  n by n matrix of quadratic coefficients
  %   B  n by 1 column of linear coefficients
  %   known  #known list of indices to known rows in Z
  %   Y  #known by cols list of fixed values corresponding to known rows in Z
  %   Optional:
  %     Aeq  meq by n list of linear equality constraint coefficients
  %     Beq  meq by 1 list of linear equality constraint constant values
  %     Aieq  mieq by n list of linear equality constraint coefficients
  %     Bieq  mieq by 1 list of linear equality constraint constant values
  %     lx n by 1 list of lower bounds [] implies -Inf
  %     ux n by 1 list of upper bounds [] implies Inf
  % Outputs:
  %   Z  n by cols solution
  %   Optional:
  %     F  struct containing all information necessary to solve a prefactored
  %     system touching only B, Y, and optionally Beq, Bieq
  %

  A = varargin{1};
  B = varargin{2};
  known = varargin{3};
  Y = varargin{4};
  Aeq = [];
  Beq = [];
  lx = [];
  ux = [];
  if nargin >= 6
    Aeq = varargin{5};
    Beq = varargin{6};
  end

  if nargin >= 8
    Aieq = varargin{7};
    Bieq = varargin{8};
  end

  if nargin >= 10
    lx = varargin{9};
    ux = varargin{10};
  end

  % number of rows
  n = size(A,1);
  % number of cols
  cols = size(Y,2);
  assert(cols == 1);

  if isempty(lx)
    lx = -Inf*ones(n,1);
  end

  LX = -speye(n,n);
  UX = speye(n,n);

  if isempty(ux)
    ux = Inf*ones(n,1);
  end

  if isempty(Aieq) && isempty(lx) && isempty(ux)
    warning('No inequality constraints found. Call min_quad_with_fixed directly');
  end

  assert(all(lx<=ux));

  old_Z = Inf*ones(n,cols);
  Z = -Inf(n,cols);
  threshold = eps;
  active_set_ieq = [];
  active_set_lx = [];
  active_set_ux = [];
  iter = 1;
  % repeat until satisfied
  while true 
    fprintf('iter: %d\n',iter);
    % keep track of last solution
    old_Z = Z;
    %% append active set constant bounds as known/fixed values
    %known_i = [known(:);active_set_lx(:);active_set_ux(:)];
    %Y_i = [Y;lx(active_set_lx);ux(active_set_ux)];
    %% append active set linear inequality constraints as *equality* constraints
    %Aeq_i = [Aeq;Aieq(active_set_ieq,:)];
    %Beq_i = [Beq;Bieq(active_set_ieq)];
    % solve equality problem
    %[Z,F,Lambda] = min_quad_with_fixed(A,B,known_i,Y_i,Aeq_i,Beq_i);
    % we need values of the lagrange multipliers for all inequality constraints
    % including bounds, so append all active set constraints as linear equality
    % constraints
    Aeq_i = [Aeq;LX(active_set_lx,:);UX(active_set_ux,:);Aieq(active_set_ieq,:)];
    Beq_i = [Beq;lx(active_set_lx)  ;ux(active_set_ux)  ; Bieq(active_set_ieq)];
    % solve equality problem
    [Z,F,Lambda] = min_quad_with_fixed(A,B,known,Y,Aeq_i,Beq_i);
    % look at lagrange multipliers of active set, remove some which are
    % negative
    Lambda_lx = Lambda(size(Aeq,1) + (1:numel(active_set_lx)));
    Lambda_ux = ... 
      Lambda(size(Aeq,1) + numel(active_set_lx) + (1:numel(active_set_ux)));
    Lambda_ieq = ... 
      Lambda( ...
        size(Aeq,1) + numel(active_set_lx) + numel(active_set_ux) + ...
        (1:numel(active_set_ieq)));
    inactive_threshold = -eps;
    active_set_lx = active_set_lx(Lambda_lx > inactive_threshold);
    active_set_ux = active_set_ux(Lambda_ux > inactive_threshold);
    active_set_ieq = active_set_ieq(Lambda_ieq > inactive_threshold);

    % Append constraints for infeasible values
    active_set_lx = [active_set_lx(:); find(Z<lx)];
    active_set_ux = [active_set_ux(:); find(Z>ux)];
    if ~isempty(Aieq)
      active_set_ieq = [active_set_ieq; find(Aieq*Z>Bieq)];
    end
    % get rid of duplicates
    active_set_lx = unique(active_set_lx);
    active_set_ux = unique(active_set_ux);
    active_set_ieq =unique(active_set_ieq);

    diff = sum((old_Z(:)-Z(:)).^2);
    fprintf('diff: %g\n',diff);
    if diff < threshold
      break;
    end 
    iter = iter + 1;
  end
end
