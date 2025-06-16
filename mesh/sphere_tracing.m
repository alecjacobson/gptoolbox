function [tI,tt,tC] = sphere_tracing(f,O,D,T,varargin)
  % [tI,tt,tC] = sphere_tracing(f,O,D,T)
  % 
  % Inputs:
  %   f  function handle f(X) : R^3 -> R signed distance function
  %   O  #O by 3 list of origin points
  %   D  #O by 3 list of direction vectors
  %   T  #O by 1 list of maximum distances
  %   Optional:
  %      Check  function handle f_chk(X) : R^3 -> R signed distance function
  %             for checking the distance function. Default is f.
  %      Lambda  inverse step size for the sphere tracing (lipschitz constant). Default is 1.
  % Outputs:
  %   tI #I list of indices into O of hit points
  %   tt #I list of parametric distances to hit points
  %   tC #I list of number of counts so that tC(i) = k, means that the k-th hit
  %     along ray tI(i) was at distance tt(i)
  %
  % Example:
  % ts = accumarray(tI(:), tt(:), [], @(x) {x});

  function [hit,t] = next_hit(f,O,D,T,t)
    for iter = 1:max_next_hit_iters
      assert(all(t<T));
      X = O + D.*(t);
      %assert(all(f(X) == f_chk(X)));
      s = abs(f(X));
      hit = s < tol;
      % don't move the hits.
      t(~hit) = t(~hit) + s(~hit)/lambda;
      escape = t > T;
      callback(struct('O',O,'D',D,'T',T,'t',t,'s',s,'hit',hit,'escape',escape));
      if any(hit | escape)
        alive = ~hit & ~escape;
        if any(alive)
          [hit(alive),t(alive)] = next_hit( ...
            f,O(alive,:),D(alive,:),T(alive),t(alive));
          return;
        else
          break;
        end
      end
    end
    if iter >= max_next_hit_iters
      max_iters_was_exceeded = true;
    end
  end

  function t = nudge(f,O,D,t)
    % they should all start as hits
    assert(all(abs(f(O + D.*t)) < tol));
    step = tol;
    for iter = 1:max_nudge_iters
      t = t + step;
      X = O + D.*t;
      %assert(all(f(X) == f_chk(X)));
      s = abs(f(X));
      unlocked = s > tol;
      if any(unlocked)
        still_hit = ~unlocked;
        if any(still_hit)
          [t(still_hit)] = nudge( ...
            f,O(still_hit,:),D(still_hit,:),t(still_hit));
        else
          break;
        end
      end
    end
  end


  function [tI,tt,tC] = recurse(f,O,D,T,C,t,I)
    assert(all(t<T));
    [hit,t] = next_hit(f,O,D,T,t);

    %hold on;
    %sct(O,'k','filled');
    %txt(O,num2str(I),'BackgroundColor',0.9*[1 1 1]);
    %qvr(O,D.*T,0,'k','ShowArrowHead','off');
    %sct(O(hit,:)+D(hit,:).*t(hit,:),'r','filled');
    %qvr(O(hit,:),D(hit,:).*t(hit,:),0,'r','ShowArrowHead','off','LineWidth',2);
    %sct(O(~hit,:)+D(~hit,:).*T(~hit,:),'b','filled');
    %hold off;
    %axis equal;


    %hold on;
    %sct(O(hit,:)+D(hit,:).*t(hit,:),'g','filled');
    %hold off;
    %pause

    if any(hit)
      C(hit) = C(hit) + 1;
      alive = hit;
      [t(alive)] = nudge(f,O(alive,:),D(alive,:),t(alive));
      alive = alive & (t < T);
      if any(alive)
        [tI,tt,tC] = recurse(f,O(alive,:),D(alive,:),T(alive),C(alive),t(alive),I(alive));
        tI = [tI;I(alive)];
        tt = [tt;t(alive)];
        tC = [tC;C(alive)];
      else
        tI = [];
        tt = [];
        tC = [];
      end
    else
      tI = [];
      tt = [];
      tC = [];
    end
  end

  % default values
  f_chk = [];
  lambda = 1;
  callback = @(local) [];
  tol = 1e-6;
  max_next_hit_iters = 100;
  max_nudge_iters = 100;
  max_iters_was_exceeded = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Check','Lambda','Callback','Tol','MaxIters'}, ...
    {'f_chk','lambda','callback','tol','max_next_hit_iters'});
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

  if isempty(f_chk)
    f_chk = f;
  end

  ts = cell(size(O,1),1);
  t = zeros(size(O,1),1);
  I = (1:size(O,1))';
  C = zeros(size(O,1),1);
  [tI,tt,tC] = recurse(f,O,D,T,C,t,I);
  tI = flip(tI);
  tt = flip(tt);
  tC = flip(tC);

  if max_iters_was_exceeded
      warning('Max iterations exceeded');
  end

end
