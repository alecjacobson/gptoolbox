function [Z,F,Lambda] = ...
  min_quad_with_fixed_active_set(varargin)
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
  %     F  see output
  % Outputs:
  %   Z  n by cols solution
  %   Optional:
  %     F  struct containing all information necessary to solve a prefactored
  %     system touching only B, Y, and optionally Beq, Bieq
  %

  %max_iter = inf;
  max_iter = 100;

  A = varargin{1};
  B = varargin{2};
  known = varargin{3};
  Y = varargin{4};
  Aeq = [];
  Beq = [];
  Aieq = [];
  Bieq = [];
  lx = [];
  ux = [];
  F = [];
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

  if nargin >= 11
    F = varargin{11};
  end

  %if isempty(F)
  %  F.Z0 = [];
  %  F.as_ieq = [];
  %  F.as_lx = [];
  %  F.as_ux = [];
  %end
  if ~isfield(F,'Z0')
    F.Z0= [];
  end
  if ~isfield(F,'as_ieq')
    F.as_ieq = [];
  end
  if ~isfield(F,'as_lx')
    F.as_lx= [];
  end
  if ~isfield(F,'as_ux')
    F.as_ux= [];
  end

  % number of rows
  n = size(A,1);
  % number of cols
  cols = size(Y,2);
  assert(cols == 1);

  if isempty(lx)
    lx = -Inf*ones(n,1);
  end

  if isempty(ux)
    ux = Inf*ones(n,1);
  end

  if isempty(Aieq) && isempty(lx) && isempty(ux)
    warning('No inequality constraints found. Call min_quad_with_fixed directly');
  end

  % Convert constant bound constraints into linear constraints before adding
  % them to the active set. Slower, but upp left corner of system matrix stays
  % the same so potentially there's room for optimization
  treat_constant_bound_as_linear_constraints = false;
  if treat_constant_bound_as_linear_constraints
    % Lower bound constraints as linear constraints, notice revert sign since
    % the <= has become a >=
    LX = -speye(n,n);
    % Upper bound constraints as linear constraints
    UX = speye(n,n);
  end

  assert(all(lx<=ux));

  old_Z = Inf*ones(n,cols);
  %Z = -Inf(n,cols);
  Z = F.Z0;
  old_Z = -Inf(n,cols);

  threshold = eps;
  %threshold = 1e-8;

  F.iter = 1;
  % repeat until satisfied
  while true 
    %fprintf('iter: %d\n',F.iter);

    if ~isempty(Z)
      % keep track of last solution
      old_Z = Z;
      % Append constraints for infeasible values
      ths = 0;
      new_as_lx = find((Z-lx)<-ths);
      new_as_ux = find((Z-ux)>ths);

      % only add k worst offenders
      min_k = 1;
      [~,ii] = sort(Z(new_as_lx),'ascend');
      k = max(min_k,ceil(numel(new_as_lx)*0.05));
      new_as_lx = new_as_lx(ii(1:min(end,k)));
      [~,ii] = sort(Z(new_as_ux),'descend');
      k = max(min_k,ceil(numel(new_as_ux)*0.05));
      new_as_ux = new_as_ux(ii(1:min(end,k)));

      F.as_lx = [F.as_lx(:); new_as_lx]; 
      F.as_ux = [F.as_ux(:); new_as_ux]; 
      if ~isempty(Aieq)
        F.as_ieq = [F.as_ieq; find(Aieq*Z>Bieq)];
      end
      % get rid of duplicates
      %F.as_lx = setdiff(unique(F.as_lx),known);
      %F.as_ux = setdiff(unique(F.as_ux),known);
      F.as_lx = unique(F.as_lx);
      F.as_ux = unique(F.as_ux);
      F.as_ieq = unique(F.as_ieq);
    end

    %fprintf('E before: %0.20f\n',energy(Z));
    if treat_constant_bound_as_linear_constraints
      % we need values of the lagrange multipliers for all inequality constraints
      % including bounds, so append all active set constraints as linear equality
      % constraints
      Aeq_i = [Aeq;LX(F.as_lx,:);UX(F.as_ux,:);Aieq(F.as_ieq,:)];
      Beq_i = [Beq;lx(F.as_lx)  ;ux(F.as_ux)  ; Bieq(F.as_ieq)];
      % solve equality problem
      [Z,mqwf,Lambda] = min_quad_with_fixed(A,B,known,Y,Aeq_i,Beq_i);
      % look at lagrange multipliers of active set, remove some which are
      % negative
      Lambda_lx = Lambda(size(Aeq,1) + (1:numel(F.as_lx)));
      Lambda_ux = ... 
        Lambda(size(Aeq,1) + numel(F.as_lx) + (1:numel(F.as_ux)));
    else
      % append active set constant bounds as known/fixed values
      known_i = [known(:);F.as_lx(:);F.as_ux(:)];
      Y_i = [Y;lx(F.as_lx);ux(F.as_ux)];
      % append active set linear inequality constraints as *equality* constraints
      Aeq_i = [Aeq;Aieq(F.as_ieq,:)];
      Beq_i = [Beq;Bieq(F.as_ieq)];
      % solve equality problem
      [Z,mqwf,Lambda,Lambda_known] = min_quad_with_fixed(A,B,known_i,Y_i,Aeq_i,Beq_i);
      % Lower bound lambda values need to be reversed because constraints
      % always represent <=
      Lambda_lx = -Lambda_known(numel(known) + (1:numel(F.as_lx)));
      Lambda_ux = ...
        Lambda_known(numel(known) + numel(F.as_lx) + ...
        (1:numel(F.as_ux)));
    end
    %fprintf('E after: %0.20f\n',energy(Z));

    Lambda_ieq = ... 
      Lambda( ...
        size(Aeq,1) + numel(F.as_lx) + numel(F.as_ux) + ...
        (1:numel(F.as_ieq)));

    %Lambda_lx
    %Lambda_ux

    %if F.iter == 2
    %  return
    %end 
    inactive_threshold = eps;


    if(isfield(F,'im'))
      AS = 0.5*ones(n,1);
      AS(F.as_lx) = 0;
      AS(F.as_ux) = 1;
      h = size(get(F.im,'CData'),1);
      w = size(get(F.im,'CData'),2);
      newAS = 0.5*ones(n,1);
      newAS(F.as_lx(Lambda_lx > inactive_threshold)) = 0;
      newAS(F.as_ux(Lambda_ux > inactive_threshold)) = 1;
      set(F.im,'CData', ...
        cat(3,reshape(AS,[h w]),reshape(newAS,[h w]),zeros([h w 1])) ...
        );
      drawnow;
      %if F.iter ==1
        %input('');
      %end
    end

    F.as_lx = F.as_lx(Lambda_lx > inactive_threshold);
    F.as_ux = F.as_ux(Lambda_ux > inactive_threshold);
    F.as_ieq = F.as_ieq(Lambda_ieq > inactive_threshold);

    diff = sum((old_Z(:)-Z(:)).^2);
    %fprintf('diff: %g\n',diff);
    if diff < threshold
      break;
    end 
    if F.iter == max_iter
      warning( ...
        sprintf('Max iterations %d reached without convergence',max_iter));
      break
    end
    F.iter = F.iter + 1;
  end

  function E = energy(ZZ)
    if isempty(ZZ)
      E = inf;
      return;
    end
    projZ = ZZ;
    projZ(ZZ>ux) = ux(ZZ>ux);
    projZ(ZZ<lx) = lx(ZZ<lx);
    if ~isempty(Aeq)
      % fix me (need to project linear inequality constraints)
      assert(false);
    end
    E = projZ'*A*projZ + projZ'*B;
  end

end
