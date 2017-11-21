function [X,Z,state] = admm(argmin_X,argmin_Z,A,B,c,state,varargin)

  % Inputs:
  %   argmin_X  Function handle returning the optimizer of:
  %      argmin  f(X) + ρ/2‖ A*X + B*Z - c + U‖²
  %        X
  %      [X,data] = argmin_X(Z,U,rho,data)
  %      Inputs:
  %        Z  #Z by dim list of dual variables
  %        U  #c by dim list of scaled Lagrange multipliers
  %        rho  current penalty parameter
  %        data  empty [] on first call, or `data` ouput from previous call
  %      Outputs:
  %        X  #X by dim list of primary variables
  %        data  persistent callback data
  %   argmin_Z  Function handle returning the optimizer of:
  %      argmin  g(Z) + ρ/2‖ A*X + B*Z - c + U‖²
  %        Z
  %      [Z,data] = argmin_Z(X,U,rho,data)
  %      Inputs:
  %        X  #X by dim list of primary variables
  %        U  #c by dim list of scaled Lagrange multipliers
  %        rho  current penalty parameter
  %        data  empty [] on first call, or `data` ouput from previous call
  %      Outputs:
  %        Z  #Z by dim list of dual variables
  %        data  persistent callback data
  %   A  #c by #X constraint matrix coefficients corresponding to rows of X
  %   B  #c by #Z constraint matrix coefficients corresponding to rows of Z
  %   c  #c by dim list of constraint constants
  %   state  (see output)
  % Outputs
  %   X  #X by dim primary part of solution
  %   Z  #Z by dim dual part of solution
  %   state  struct containing persistent/reusable data
  %  
  %   

  max_iter = 2000;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MaxIter'}, ...
    {'max_iter'});
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

  tol_abs = 1e-8;
  tol_rel = 1e-6;
  check_interval = 10;
  bmu = 5;
  btao_inc = 2;
  btao_dec = 2;

  % Initial conditions
  if isempty(state)
    state.X = rand(size(A,2),size(c,2));
    state.Z = rand(size(B,2),size(c,2));
    state.U = zeros(size(c,1),size(c,2));
    state.rho_prev = nan;
    state.rho = 1;
    state.argmin_X_data = [];
    state.argmin_Z_data = [];
  end

  cnorm = norm(c,'fro');
  for iter = 1:max_iter
    [state.X,state.argmin_X_data] = argmin_X(state.Z,state.U,state.rho,state.argmin_X_data);
    state.Z_prev = state.Z;
    [state.Z,state.argmin_Z_data] = argmin_Z(state.X,state.U,state.rho,state.argmin_Z_data);
    state.U_prev = state.U;
    state.U = state.U+A*state.X+B*state.Z-c;
    state.rho_prev = state.rho;
      %Sprev = state.Z_prev(1:4780,:);
      %Eprev = state.Z_prev(4780+(1:4780),:);
      %S = state.Z(1:4780,:);
      %E = state.Z(4780+(1:4780),:);
      %dual_residual = sqrt(state.rho*sum((S(:)-Sprev(:)).^2+(E(:)-Eprev(:)).^2))
    dual_residual = state.rho*norm(A'*B*(state.Z_prev - state.Z),'fro');
    residual = norm(A*state.X+B*state.Z-c,'fro');
    if mod(iter,check_interval) == 0
      if residual > bmu*dual_residual
        state.rho = btao_inc*state.rho;
        % From python code: https://github.com/tneumann/cmm/blob/master/cmmlib/cmm.py
        state.U = state.U/btao_inc;
      elseif dual_residual > bmu*residual
        state.rho = state.rho/btao_dec;
        state.U = state.U*btao_dec;
      end
    end
    % From Python code
    k = size(c,2);
    eps_pri = sqrt(k*2)*tol_abs + tol_rel*max([norm(A*state.X,'fro'),norm(B*state.Z,'fro'),cnorm]);
    eps_dual = sqrt(k)*tol_abs +  tol_rel*state.rho*norm(A'*state.U,'fro');
    if residual < eps_pri && dual_residual < eps_dual
      break;
    end
  end
  X = state.X;
  Z = state.Z;
end
