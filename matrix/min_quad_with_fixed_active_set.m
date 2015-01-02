function [Z,F,Lambda] = min_quad_with_fixed_active_set( ...
    A,B,known,Y,Aeq,Beq,Aieq,Bieq,lx,ux,F,varargin)
  % MIN_QUAD_WITH_FIXED_ACTIVE_SET Minimize quadratic energy Z'*A*Z + Z'*B + C
  % with constraints that Z(known) = Y, optionally also subject to the
  % constraints Aeq*Z = Beq, and further optionally subject to the linear
  % inequality constraints that Aieq*Z <= Bieq and constant inequality
  % constraints lx <= x <= ux
  %
  % [Z] = min_quad_with_fixed_active_set(A,B,known,Y,Aeq,Beq,Aeiq,Beiq,lx,ux)
  % [Z,F,Lambda] = ...
  %   min_quad_with_fixed_active_set(A,B,known,Y,Aeq,Beq,Aeiq,Beiq,lx,ux,F)
  %
  % Inputs:
  %   A  n by n matrix of quadratic coefficients
  %   B  n by 1 column of linear coefficients
  %   known  #known list of indices to known rows in Z
  %   Y  #known by cols list of fixed values corresponding to known rows in Z
  %   Optional:
  %     Aeq  meq by n list of linear equality constraint coefficients
  %     Beq  meq by 1 list of linear equality constraint constant values
  %     Aieq  mieq by n list of linear inequality constraint coefficients
  %     Bieq  mieq by 1 list of linear inequality constraint constant values
  %     lx n by 1 list of lower bounds [] implies -Inf
  %     ux n by 1 list of upper bounds [] implies Inf
  %     F  see output
  %     'Quiet' followed by whether to suppress output
  %     'MaxIter' followed by maximum number of iterations {100}
  %     'ActiveThreshold' followed by threshold to determine whether constraint
  %       needs to be activated (smaller and negative means more
  %       agressive, larger means more lenient) {+eps}
  %     'InactiveThreshold' followed by threshold to determine whether
  %       constraint should stay active. Smaller and negative means constraints
  %       are "stickier", larger means constraints flicker more {+eps}
  %     'StopThreshold' followed by threshodl to determine whether no more
  %       progress is being made based on difference between iterations'
  %       solutions.
  % Outputs:
  %   Z  n by cols solution
  %   Optional:
  %     F  struct containing all information necessary to solve a prefactored
  %     system touching only B, Y, and optionally Beq, Bieq
  %
  % Examples:
  %   qZ = quadprog( ... 
  %     2*A,B,Aieq,Bieq, ...
  %     [Aeq;sparse(1:numel(known),known,1,numel(known),size(A,2))], ...
  %     [Beq;Y],lx,ux);
  %
  %   % bbw works well with 'InactiveThreshold',-1e-9,'ActiveThreshold',0
  %

  %A = varargin{1};
  %B = varargin{2};
  %known = varargin{3};
  %Y = varargin{4};
  %Aeq = [];
  %Beq = [];
  %Aieq = [];
  %Bieq = [];
  %lx = [];
  %ux = [];
  %F = [];
  %if nargin >= 6
  %  Aeq = varargin{5};
  %  Beq = varargin{6};
  %end

  %if nargin >= 8
  %  Aieq = varargin{7};
  %  Bieq = varargin{8};
  %end

  %if nargin >= 10
  %  lx = varargin{9};
  %  ux = varargin{10};
  %end

  %if nargin >= 11
  %  F = varargin{11};
  %end

  %if isempty(F)
  %  F.Z0 = [];
  %  F.as_ieq = [];
  %  F.as_lx = [];
  %  F.as_ux = [];
  %end
  if ~exist('F','var')
      F = [];
  end
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
  max_iter = 100;
  inactive_threshold = eps;
  stop_threshold = eps;
  active_threshold = eps;
  quiet = false;

  % solve equality problem
  force_Aeq_li = [];

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Quiet','MaxIter','ForceAeqLI','InactiveThreshold','StopThreshold','ActiveThreshold'}, ...
    {'quiet','max_iter','force_Aeq_li','inactive_threshold','stop_threshold','active_threshold'});
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


  % number of rows
  n = size(A,1);
  % number of cols
  if isempty(Y)
    cols = 1;
  else
    cols = size(Y,2);
    assert(cols == 1);
  end
  
  if numel(lx) == 1
    lx = lx*ones(n,1);
  end
  
  if numel(ux) == 1
    ux = ux*ones(n,1);
  end

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

  F.iter = 1;
  % repeat until satisfied
  while true 
    %fprintf('iter: %d\n',F.iter);
    %tic;

    if ~isempty(Z)
      % Append constraints for infeasible values
      new_as_lx = find((Z-lx)<-active_threshold);
      new_as_ux = find((Z-ux)>active_threshold);
      new_as_ieq = [];
      if ~isempty(Aieq)
        new_as_ieq = find((Aieq*Z-Bieq)>active_threshold);
      end

      %fprintf('%0.17f %d %d %g\n', ...
      %  energy(Z),numel(new_as_lx), numel(new_as_ux),max(max(lx-Z,Z-ux)));
      diff = sum((old_Z(:)-Z(:)).^2);
      %fprintf('diff: %g\n',diff);
      if diff < stop_threshold
        if numel(new_as_lx) > 0 || numel(new_as_ux) > 0
          % Indeed this is "disturbing" at first glance. Why should there be
          % new constraints? Turns out these "last constraints" are only
          % producing numerical noise. The can safely be caught and ignored by
          % setting active_threshold to eps rather than 0
          %warning( ...
          %  sprintf('Diff (%g) < threshold (%g) but new constraints (%d)', ...
          %    diff,threshold,numel(new_as_lx)+numel(new_as_ux)));
        end
        break;
      end 
      % keep track of last solution
      old_Z = Z;

      %% only add k worst offenders
      %min_k = 1;
      %[~,ii] = sort(Z(new_as_lx),'ascend');
      %k = max(min_k,ceil(numel(new_as_lx)*0.05));
      %new_as_lx = new_as_lx(ii(1:min(end,k)));
      %[~,ii] = sort(Z(new_as_ux),'descend');
      %k = max(min_k,ceil(numel(new_as_ux)*0.05));
      %new_as_ux = new_as_ux(ii(1:min(end,k)));

      F.as_lx = [F.as_lx(:); new_as_lx]; 
      F.as_ux = [F.as_ux(:); new_as_ux]; 
      F.as_ieq = [F.as_ieq(:); new_as_ieq];
      % get rid of duplicates
      %F.as_lx = setdiff(unique(F.as_lx),known);
      %F.as_ux = setdiff(unique(F.as_ux),known);
      F.as_lx = unique(F.as_lx);
      F.as_ux = unique(F.as_ux);
      F.as_ieq = unique(F.as_ieq);
      F.as_lx = setdiff(F.as_lx,known);
      F.as_ux = setdiff(F.as_ux,known);
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
      Lambda_ieq = ... 
        Lambda( ...
          size(Aeq,1) + numel(F.as_lx) + numel(F.as_ux) + ...
          (1:numel(F.as_ieq)));
    else
      % append active set constant bounds as known/fixed values
      known_i = [known(:);F.as_lx(:);F.as_ux(:)];
      Y_i = [Y;lx(F.as_lx);ux(F.as_ux)];
      % append active set linear inequality constraints as *equality* constraints
      Aeq_i = [Aeq;Aieq(F.as_ieq,:)];
      Beq_i = [Beq;Bieq(F.as_ieq)];
      [Z,~,Lambda,Lambda_known] = min_quad_with_fixed(A,B,known_i,Y_i,Aeq_i,Beq_i,struct('force_Aeq_li',force_Aeq_li));
      %[Z,~,Lambda,Lambda_known] = min_quad_with_fixed(A,B,known_i,Y_i,Aeq_i,Beq_i,struct('force_Aeq_li',true));
      %[Z_qr,~,Lambda_qr,Lambda_known] = min_quad_with_fixed(A,B,known_i,Y_i,Aeq_i,Beq_i,struct('force_Aeq_li',false));
      %save('bad.mat','A','B','known_i','Y_i','Aeq_i','Beq_i');
      %fprintf('Z: %g Lambda: %g\n', ...
      %  max(max(abs(Z-Z_qr))), ...
      %  max(max(abs(Lambda-Lambda_qr))));
      %  input('');

      % Lower bound lambda values need to be reversed because constraints
      % always represent <=
      Lambda_lx = -Lambda_known(numel(known) + (1:numel(F.as_lx)));
      Lambda_ux = ...
        Lambda_known(numel(known) + numel(F.as_lx) + ...
        (1:numel(F.as_ux)));
      Lambda_ieq = Lambda(size(Aeq,1) + (1:numel(F.as_ieq)));
    end
    %fprintf('E after: %0.20f\n',energy(Z));

    %Lambda_lx
    %Lambda_ux

    %if F.iter == 2
    %  return
    %end 


    if isfield(F,'im') || isfield(F,'t')
      AS = 0.5*ones(n,1);
      AS(F.as_lx) = 0;
      AS(F.as_ux) = 1;
      newAS = 0.5*ones(n,1);
      newAS(F.as_lx(Lambda_lx > inactive_threshold)) = 0;
      newAS(F.as_ux(Lambda_ux > inactive_threshold)) = 1;
      if isfield(F,'t')
        set(F.t,'FaceVertexCData',[AS newAS zeros(n,1)]);
        %set(F.t,'CData',log10(abs(Z)));
        colorbar;
      end
      if isfield(F,'im')
        h = size(get(F.im,'CData'),1);
        w = size(get(F.im,'CData'),2);
        set(F.im,'CData', ...
          cat(3,reshape(AS,[h w]),reshape(newAS,[h w]),zeros([h w 1])) ...
          );
      end
      drawnow;
      %if F.iter ==1
        %input('');
      %end
    end


    F.as_lx = F.as_lx(Lambda_lx > inactive_threshold);
    F.as_ux = F.as_ux(Lambda_ux > inactive_threshold);
    F.as_ieq = F.as_ieq(Lambda_ieq > inactive_threshold);

    if F.iter == max_iter
      if ~quiet
        warning( ...
          sprintf('Max iterations %d reached without convergence',max_iter));
      end
      break
    end
    F.iter = F.iter + 1;
    %toc;
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
