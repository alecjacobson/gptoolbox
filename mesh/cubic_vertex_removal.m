function [C,t,err] = cubic_vertex_removal(C1,C2,varargin)
  % CUBIC_VERTEX_REMOVAL Given a G¹ continuous sequence of cubic Bézier curves,
  % optimize the positions of a new single curve to minimize its integrated
  % squared distance to those curves. This requires both estimating the
  % parameterization mapping between the inputs and output (non-linear problem)
  % and then computing the optimal positions as a linear least squares problem.
  % 
  % C = cubic_vertex_removal(C1,C2)
  % [C,t,err] = cubic_vertex_removal(C1,C2,'ParameterName',ParameterValue,…)
  %
  % Inputs:
  %   C1  4 by dim list of first curves control point positions
  %   C2  4 by dim list of second curves control point positions, it's assumed
  %     that C1(4,:) == C2(1,:) and C1(4,:) - C1(3,:) = s * (C2(2,:) - C2(1,:))
  %     for some s>= 0.
  %     Optional:
  %       'Method' followed by one of the following:
  %         {'perfect'}  use root finding. This is fastest when you only care about
  %           finding a perfect fit. It may lead to a very poor approximate when
  %           a perfect fit is not possible. A perfect fit is when 
  %           `[C1,C2] = cubic_split(C,t)`. If this is the case, then the
  %           returned `err` for this method should be very close to zero. (Note
  %           that "perfect fit" → err=0 but err=0 does not necessarily imply
  %           "perfect fit". 
  %         'iterative'  use iterative method. This is slower but more accurate
  %           when a perfect fit is not possible.
  %       'MaxIter'  followed by maximum number of iterations for iterative
  %         method {100}
  %       't0'  followed by initial guess of t for iterative method {relative
  %         approximate arc lengths}
  %       'AlreadyGenerated'  followed by whether the automatically generated helper
  %         functions cubic_vertex_removal_polyfun and cubic_vertex_removal_g
  %         have already been generated. On some machines checking `exist()` can
  %         be really slow. So, you could call
  %         `cubic_vertex_removal(…,'AlreadyGenerated',false)` onces to
  %         generate the files and the call
  %         `cubic_vertex_removal(…,'AlreadyGenerated',true) for subsequent
  %         calls {false}.
  % Outputs:
  %   C  4 by dim list of output coordinates. By default:
  %     C(1,:) = C1(1,:),
  %     C(4,:) = C2(4,:), (C₀ continuity)
  %       and
  %     C(2,:) - C(1,:) = s1*(C1(2,:) - C1(1,:)) with s1>=0
  %     C(4,:) - C(3,:) = s2*(C2(4,:) - C2(3,:)) with s2>=0 (G₁ continutity)
  %   t scalar value between [0,1] defining the piecewise-linear
  %     parameterization mapping between the inputs and output.  C1 is mapped to
  %     [0,t] and C2 is mapped to [t1,1].
  %   err  Integrated squared distance between the output and input curves (see
  %   cubic_cubic_integrated_distance).
  %
  % Example:
  %   C = [0 0;1 1;2 -1;3 0];
  %   tgt = 0.1;
  %   [C1,C2] = cubic_split(C,tgt);
  %   [C,t,err] = cubic_vertex_removal(C1,C2,'Method','perfect');
  %   clf;
  %   hold on;
  %   plot_cubic(C1,[],[],'Color',orange);
  %   plot_cubic(C2,[],[],'Color',orange);
  %   plot_cubic(C,[],[],'Color',blue);
  %   hold off;
  %   axis equal;
  %   set(gca,'YDir','reverse')
  %   title(sprintf('err: %g',err),'FontSize',30);
  %   

  method = 'perfect';
  max_iter = 100;
  t0 = [];
  promise_already_built = false;
  E_tol = 1e-15;
  grad_tol = 1e-8;

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Method','MaxIter','t0','AlreadyGenerated','Tol','GradientTol',}, ...
    {'method','max_iter','t0','promise_already_built','E_tol','grad_tol'});
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

  if strcmp(method,'cubic-polish')
    [C,t,err] = cubic_vertex_removal(C1,C2,varargin{:},'Method','cubic');
    % Slip in 'MaxIter',1 before user parameters so it gets over-written if
    % provided. 
    [C,t,err] = cubic_vertex_removal(C1,C2,'MaxIter',1,varargin{:},'t0',t,'Method','iterative');
    return;
  end


  switch method
  case {'perfect','cubic'}
    % _If_ t is perfect, then C is a simple function of t.
    C_from_t = @(C1,C2,t) ...
      [C1(1,:); ...
          (1./t)*(C1(2,:)-C1(1,:)) + C1(1,:); ...
      (1./(1-t))*(C2(end-1,:)-C2(end,:)) + C2(end,:); ...
      C2(end,:)];
    switch method
    case 'perfect'
      % Build a root finder for the 1D problem
      if ~promise_already_built && ~exist('cubic_vertex_removal_polyfun','file');
        warning('assuming L2 not l2 error');
        dim = 1;
        % Matlab is (sometimes?) confused that this is a static workspace and
        % refuses to let syms create variables.
        %syms('iC1',[4 dim],'real');
        %syms('iC2',[4 dim],'real');
        %% complains if I mark st as 'real'
        %syms('st',[1 1]);
        iC1 = [
          sym('iC11','real')
          sym('iC12','real')
          sym('iC13','real')
          sym('iC14','real')];
        iC2 = [
          sym('iC21','real')
          sym('iC22','real')
          sym('iC23','real')
          sym('iC24','real')];
        st = sym('st');

        sC = C_from_t(iC1,iC2,st);
        [oC1,oC2] = cubic_split(sC,st);
        % L2 style.
        sg = sum( [oC1-iC1;oC2-iC2].^2, 'all');
        sres = solve(diff(sg,st) == 0,st);
        %sdgdt = simplify(diff(sg,st)); 
        %vroots = @(C1,C2) double(vpasolve(subs(subs(sdgdt,iC1,C1),iC2,C2)==0,st,[0 1]));
        sres_children = children(sres(1));
        % Should cache these two:
        %polyfun = matlabFunction(flip(coeffs(sres_children(1),sres_children(2))),'Vars',{iC1,iC2});
        % Matlab2023a seems to use slightly different cell arrays for sres_children. 
        polyfun = matlabFunction( ...
          flip(coeffs(sres_children{1},sres_children{2})), ...
          'Vars',{iC1,iC2}, ...
          'File','cubic_vertex_removal_polyfun');
        g = matlabFunction(sg,'Vars',{iC1,iC2,st},'File','cubic_vertex_removal_g');
      else
        polyfun = @cubic_vertex_removal_polyfun;
        g = @cubic_vertex_removal_g;
      end
      keepreal = @(C) C(imag(C)==0 & real(C)>=0 & real(C)<=1);
      % Using numerical roots is faster than vroots
      nroots = @(K1,K2) keepreal(roots(polyfun(K1,K2)));
      keepmin = @(ts,Es) ts(find(Es==min(Es),1));
      % We'll determine t based on the 1D g function. But we'll compute the
      % returned energy below using the full dim-D problem.
      keepmin_g = @(K1,K2,ts) keepmin(ts,arrayfun(@(t) g(K1,K2,t),ts));
      find_t = @(K1,K2) keepmin_g(K1,K2,nroots(K1,K2));
      % Decide which coordinate to use (pick a non-degenerate one).
      % Based on max-extent.
      [~,i] = max(max([C1(1,:);C2(end,:)])-min([C1(1,:);C2(end,:)]));
      t = find_t(C1(:,i),C2(:,i));
      %assert(~isempty(t));
      if isempty(t)
        C = nan(4,2);
        t = nan;
        err = inf;
        return;
      end
    case 'cubic'
      D1 = -6.*C1(1,:) + 18.*C1(2,:) - 18.*C1(3,:) + 6.*C1(4,:);
      D2 = -6.*C2(1,:) + 18.*C2(2,:) - 18.*C2(3,:) + 6.*C2(4,:);
      D2_sqr_len = sum(D2.*D2,2);
      D1_sqr_len = sum(D1.*D1,2);
      D2_D1 = sum(D2.*D1,2);
      r = D2_D1/D1_sqr_len;
      t = (1 + (r.^(1/3))).^-1;
    end
    C = C_from_t(C1,C2,t);
  case 'iterative'
    % use arc-length to guess t₁
    if isempty(t0)
      tol = 1e-5;
      ts = matrixnormalize(cumsum([0;spline_arc_lengths([C1;C2],[1 2 3 4;5 6 7 8],tol)]));
      t0 = ts(2);
    end
    t1 = t0;
    % Build null space matrices. so that C(:) = S*V + B(:) satisfies C₀ and G₁
    % constraints for any V
    B = [C1(1,:);C1(1,:);C2(4,:);C2(4,:)];
    B1 = [0 0;(C1(2,:) - C1(1,:));0 0;0 0];
    B2 = [0 0;0 0;(C2(3,:) - C2(4,:));0 0];
    S = [B1(:) B2(:)];

    f = @(t1) objective_t1(C1,C2,t1,B,S);
    [E,C] = f(t1);
    for iter = 1:max_iter
      %dfdt1 = (f(t1+1e-5)-f(t1-1e-5))/(2*1e-5);
      % Complex step is bit faster and more accurate/stable. For 1D input, I
      % believe this should be as good as autodiff and probably as good as we
      % can get without a lot of hand derivativation/compiling code.
      dfdt1 = imag(f(complex(t1,1e-100)))/1e-100;
      
      if norm(dfdt1,inf) < grad_tol
        break;
      end
      dt1 = 0.5*sign(-dfdt1);
      [alpha,t1] = backtracking_line_search(f,t1,dfdt1,dt1,0.3,0.5);
      if alpha == 0
        %warning('line search failed');
        break;
      end
      [E,C] = f(t1);
      if E < E_tol
        break;
      end
    end

    t = t1;
  otherwise
    error(['unknown method :' method]);
  end


  if nargout>2
    % L2 not l2
    [E1] = cubic_cubic_integrated_distance( ...
      0,t, ...
      1/t,0, ...
      C1, ...
      1,0, ...
      C);
    [E2] = cubic_cubic_integrated_distance( ...
      t,1, ...
      1/(1-t),-t/(1-t), ...
      C2, ...
      1,0, ...
      C);
    err = E1+E2;
  end


  % Helper functions for iterative method
  function [E,C] = objective_t1(C1,C2,t1,B,S)
    if isfloat(t1) && (t1>1 || t1<0)
      E = inf;
      return;
    end
    C = nan(4,2);
    [~,H,F,c] = objective(C1,C2,C,t1);
    % Enforce constraints via subspace
    % C = [
    %   C1(1,:)
    %   C1(1,:) + v1 * (C1(2,:) - C1(1,:))
    %   C2(4,:) + v2 * (C2(3,:) - C2(4,:))
    %   C2(4,:)
    %   ];
    HH = repdiag(H,2);
    V = ((S.'*HH*S)\(-S.'*F(:)-S.'*HH*B(:)));
    V = max(V,0);
    C = reshape( B(:) + S*V , size(C));
    [E] = objective(C1,C2,C,t1);
  end
  function [E,H,F,c] = objective(C1,C2,C,t1)
    if isfloat(t1) && (t1>1 || t1<0)
      E = inf;
      return;
    end
    % WARNING
    w1 = 1;
    w2 = 1;
    %w1 = t1; 
    %w2 = 1-t1;
    if nargout == 4
      % Given t1 update C
      [H1,F1,c1,E1] = cubic_cubic_integrated_distance( ...
        0,t1, ...
        1/t1,0, ...
        C1, ...
        1,0, ...
        C);
      [H2,F2,c2,E2] = cubic_cubic_integrated_distance( ...
        t1,1, ...
        1/(1-t1),-t1/(1-t1), ...
        C2, ...
        1,0, ...
        C);
      % Equal weighting ("L2")
      H = w1*H1 + w2*H2;
      F = w1*F1 + w2*F2;
      c = w1*c1 + w2*c2;
    else 
      % Given t1 update C
      [E1] = cubic_cubic_integrated_distance( ...
        0,t1, ...
        1/t1,0, ...
        C1, ...
        1,0, ...
        C);
      [E2] = cubic_cubic_integrated_distance( ...
        t1,1, ...
        1/(1-t1),-t1/(1-t1), ...
        C2, ...
        1,0, ...
        C);
    end
    E = w1*E1 + w2*E2;
  end
end
