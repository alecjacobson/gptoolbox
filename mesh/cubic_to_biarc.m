function [G,O1,O2,R1,R2,ccw1,ccw2,valid,sol] = cubic_to_biarc(C,varargin)


  method = 'optimal';
  search_method = 'golden';
  search_method = 'uniform';
  budget = 64;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Method','SearchMethod'}, ...
    {'method','search_method'});
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

  P1 = C(1,:);
  C1 = C(2,:);
  P2 = C(4,:);
  C2 = C(3,:);

  % un-normalized
  [T1,T2] = spline_robust_endpoint_tangent_directions(C,[1 2 3 4]);

  switch method
  case {'juckett','optimal'}
    p1 = P1';
    p2 = P2';
    t1 = normalizerow(T1)';
    t2 = normalizerow(T2)';
    sol = biarc_solve_2d(p1,p2,t1,t2,false,0);
    %switch method
    %case {'juckett'}
    %  sol
    %end
  case 'incenter'
    nc = 1;
    [Q,valid,t1,t2] = intersect_lines_2d(P1,T1,P2,-T2);
    G = incenter([P1;P2;Q],reshape(1:3*nc,[],3));
    [O1,valid1,R1,~,ccw1] = circle_center_P1G_tangent(P1,P1+T1,G);
    [O2,valid2,R2,~,ccw2] = circle_center_P1G_tangent(P2,P2 - T2,G);
    %if t1 < 0
    %  ccw1 = ~ccw1;
    %end
    %if t2 < 0
    %  ccw2 = ~ccw2;
    %end
    % flip because we did (P2,C2)→G instead of G→(P2,C2)
    ccw2 = ~ccw2;
  end

  nsamples = 10;
  switch method
  case 'optimal'
    f = @(d1) biarc_cubic_mse(p1,p2,t1,t2,d1,C,nsamples);
    def_sol = sol;
    best_sol = def_sol;
    best_error = f(def_sol.d1);
    I = biarc_d1_intervals(p1,p2,t1,t2);
    %I = [0.56 1.9]*def_sol.d1;
    I = [0.125 8.0]*def_sol.d1;
    %if isinf(I(2))
    %  I(2) = 4*def_sol.d1;
    %end

    % number of function evals
    switch search_method
    case 'uniform'
      for d1 = linspace(I(1),I(2),budget)
        [err,sol] = f(d1);
        if err < best_error
          best_error = err;
          best_sol = sol;
        end
      end
    case 'golden'
      a = I(1);
      b = I(2);
      phi = (1 + sqrt(5)) / 2; % golden ratio
      c = b - (b - a) / phi;
      d = a + (b - a) / phi;
      for i = 1:ceil(budget/2)
        [err_c, sol_c] = f(c);
        [err_d, sol_d] = f(d);
        if err_c < err_d
          b = d;
          d = c;
          c = b - (b - a) / phi;
          if err_c < best_error
            best_error = err_c;
            best_sol = sol_c;
          end
        else
          a = c;
          c = d;
          d = a + (b - a) / phi;
          if err_d < best_error
            best_error = err_d;
            best_sol = sol_d;
          end
        end
      end
    end
    %fprintf('%g $%g | %g %g %g: $%g\n',def_sol.d1, biarc_cubic_mse(p1,p2,t1,t2,def_sol.d1,C,nsamples), ...
    %I(1),best_sol.d1,I(2),best_error);
    %
    sol = best_sol;
  end

  switch method
  case {'juckett','optimal'}
    G = sol.pm';
    O1 = sol.c1';
    O2 = sol.c2';
    R1 = sol.r1;
    R2 = sol.r2;
    valid = true;
    % These can be wrong
    ccw1 = sol.theta1>0;
    ccw2 = sol.theta2>0;
  end

end
