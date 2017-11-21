function [X,state] = quadprog_box(H,f,Aeq,Beq,lx,ux,state,varargin)
  % QUADPROG_BOX Minimize a box constrained quadratic program using ADMM.
  % 
  %  min ½ trace(X'*H*X) + trace(X'*f) 
  %  subject to Aeq * X = Beq
  %         and lx ≤ X ≤ ux
  %
  % [X,state] = quadprog_box(H,f,Aieq,Bieq,Aeq,Beq,lx,ux,X0,state)
  % 
  % Inputs:
  %   H  #X by #X sparse quadratic coefficients matrix (with іmplied ½, like
  %     quadprog, but unlike min_quad_with_fixed)
  %   f  #X by dim linear quadratic coefficients vector/matrix
  %%   Aieq  ignored (todo/match quadprog api)
  %%   Bieq  ignored (todo/match quadprog api)
  %   Aeq  #Aeq by #X sparse linear equality constraints matrix
  %   Beq  #Beq by dim linear equality constraints right hand side
  %   lx  #X by dim lower bounds
  %   ux  #X by dim upper bounds
  %   state  precomputation struct (see output) {[]}
  %     X0  #X by dim initial guess {[]}
  % Outputs:
  %   X  #X by dim solution
  %   state  struct containing precomputation information. This is reusable iff
  %     the H and Aeq do not change. It is OK, for f, Beq, lx and ux to change.
  % 
  % 

  % Todo: what if H is constant but Aeq is changing? more dual variables for
  % Aeq*X = beq? For sparse_eigs(...,'brandt') this would mean a faster
  % quadratic solve (though who knows about overall convergence).
  
  max_iter = 800;
  Aeq_li = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MaxIter','Aeq_li'}, ...
    {'max_iter','Aeq_li'});
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

  function [X,data] = argmin_X(Z,U,rho,data)
    if isempty(data)
      data.rho = inf;
    end
    data.rho_prev = data.rho;
    data.rho = rho;
    if data.rho ~= data.rho_prev
      if isempty(Aeq)
        A = data.rho/2*I + H;
        [cL,p,cS] = chol(A,'lower');
        cU = cL';
        cP = cS';
        cQ = cS;
        data.X_solve = @(rhs,Beq) cQ * (cU \ (cL \ ( cP * rhs)));
      else
        mqwf = [];
        if ~isempty(Aeq_li)
          mqwf.force_Aeq_li = Aeq_li;
        end
        [~,mqwf] = min_quad_with_fixed(data.rho/2*I + 2*H,-f,[],[],Aeq,Beq,mqwf);
        data.X_solve = @(rhs,Beq) min_quad_with_fixed([],-2*rhs,[],[],[],Beq,mqwf);
      end
    end
    X = data.X_solve(-f+rho/2*((Z-U)),Beq);
  end

  function [Z,data] = argmin_Z(X,U,rho,data)
    Z = min(max(U+X,lx),ux);
  end

  
  n = size(H,1);
  I = speye(n);
  k = max([size(f,2) size(lx,2) size(ux,2) size(Beq,2)]);
  [X,~,state] = admm(@argmin_X,@argmin_Z,I,-I,zeros(n,k),state,'MaxIter',max_iter);


end
