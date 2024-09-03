function [T,F,U] = cubic_uniformly_sample(C,n,varargin)
  % [T,L,U] = cubic_uniformly_sample(C,n)
  %
  % Inputs:
  %  C  #C by 4 list of cubic control points
  %  n  number of samples
  %    Optional:
  %      'NumQuadrature' followed by number of quadrature points {10}
  %      'MaxIter' followed by maximum number of iterations {10}
  %      'Tol' followed by tolerance for function value {1e-16}
  %      'GradTol' followed by tolerance for gradient {1e-16}
  % Outputs:
  %  T  n list of parameter values
  %  L  n-1 list of segment lengths
  %  U  n by size(C,2) list of uniformly sampled points
  %
  % Example:
  %  % padded samples so that each *sample* corresponds to *center* of an equal
  %  % length segment.
  %  % (perhaps consider just using [½ 1 1 … 1 1 ½] weights instead)
  %  [T,L,P] = cubic_uniformly_sample(C,2*n+1);
  %  T = T(2:2:end);
  %  L = L(1:2:end) + L(2:2:end);
  %  P = P(2:2:end,:);

  % number of quadrature points
  f_tol = 1e-16;
  grad_tol = 1e-16;
  nq = 10;
  max_iter = 10;


  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'NumQuadrature','MaxIter','Tol','GradTol'}, ...
    {'nq','max_iter','f_tol','grad_tol'});
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

  t = linspace(0,1,n)';
  JI = repmat((1:n-1),2,1)';
  JJ = [1:n-1;2:n]';
  for iter = 1:max_iter
    [F,dFda,dFdb] = cubic_arc_length(C,nq,t(1:end-1),t(2:end));
    J = sparse(JI,JJ,[dFda dFdb],n-1,n);
    h = mean(F);
    dfdt = J'*(F-h);
    f = @(t) 0.5*sum((cubic_arc_length(C,nq,t(1:end-1),t(2:end)) - h).^2,'all');

    if norm(dfdt,inf) < grad_tol
      warning('converged');
      break;
    end

    method = 'gd';
    method = 'sobolev';
    method = 'gn';
    switch method
    case 'gd'
      dt = -dfdt;
    case 'sobolev'
      dt([1 end]) = 0;
      E = [1:numel(t)-1;2:numel(t)]';
      A = adjacency_matrix(E);
      L = diag(sum(A,2))-A;
      dt = [0;-(L(2:end-1,2:end-1) \ dfdt(2:end-1));0];
    case 'gn'
      H = J'*J;
      dt = [0;-(H(2:end-1,2:end-1) \ dfdt(2:end-1));0];
    end

    [s,t,ft] = backtracking_line_search(f,t,dfdt,dt,0.3,0.5);

    %clf;
    %plot_cubic(C,[],[],{{'LineWidth',1},{'LineWidth',0.5}});
    %hold on;
    %sct(cubic_eval(C,t),'or','filled');
    %hold off;
    %axis equal;
    %fprintf('iter %d: %g %g\n',iter,s,f(t));
    %drawnow;

    if ft<f_tol
      % recompute after step one last time
      F = cubic_arc_length(C,nq,t(1:end-1),t(2:end));
      break;
    end
    if s == 0
      error('line search failed');
    end
  end
  T = t;

  if nargout>2
    U = cubic_eval(C,t);
  end

end
