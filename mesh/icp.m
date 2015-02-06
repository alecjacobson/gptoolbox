function [R,t,BRt,e,KDTA,KDTB] = icp(A,B,varargin)
  % ICP Find a rigid transformation (R,t) that registers point set B to point
  % set A using iterative closest point (icp).
  %
  % [R,t] = icp(A,B,varargin)
  % [R,t,BRt,e,KDTA,KDTB] = icp(A,B,varargin)
  %
  % Inputs:
  %   A  #A by dim list of points
  %   B  #B by dim list of points
  %   Optional:
  %     'MaxIter' followed by maximum number of iterations {inf}
  %     'MaxRNorm'  followed by minimum norm on diff entries of R {1e-5}
  %     'MaxTNorm'  followed by minimum norm on diff entries of t {1e-5}
  %     'MaxSamples'  followed by number of random samples to use from each
  %       point set {inf --> #A,#B}
  %     'KDTreeA'  followed by kdtree on A (see output)
  %     'KDTreeB'  followed by kdtree on B (see output)
  %     'R0'  followed by dim by dim initial rotation matrix
  %     't0'  followed by 1 by dim initial rotation matrix
  %     'Quiet'  followed by whether to be quite {false}
  % Outputs:
  %   R  dim by dim rotation matrix
  %   t  1 by dim translation
  %   BRt  #B by dim list of deformed points in B, BRt = bsxfun(@plus,B*R,t)
  %   e  ICP energy
  %   KDTA  kdtree around A (see knnsearch)
  %   KDTB  kdtree around B
  %   

  % default values
  max_iter = 100;
  min_R_norm = 1e-5;
  min_t_norm = 1e-5;
  max_samples = inf;
  KDTA = [];
  KDTB = [];
  tsurf_t = [];
  R = [];
  t = [];
  quiet = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MaxIter','MinRNorm','MinTNorm','MaxSamples','tsurfHandle','KDTreeA', ...
    'KDTreeB','R0','t0','Quiet' ...
    }, ...
    {'max_iter','min_R_norm','min_t_norm','max_samples','tsurf_t','KDTA', ...
    'KDTB','R','t','quiet' ...
    });
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

  max_samples = min([max_samples;size(A,1);size(B,1)]);

  dim = size(A,2);
  bbd = norm(max([A;B])-min([A;B]));

  if isempty(R)
    R = eye(dim);
  end
  if isempty(t)
    t = zeros(1,dim);
  end
  BRt = bsxfun(@plus,B*R,t);
  % Pre-build KDtree search data structures
  if isempty(KDTA)
    KDTA = createns(A,'nsmethod','kdtree');
  end
  if isempty(KDTB)
    KDTB = createns(B,'nsmethod','kdtree');
  end

  IA = 1:size(A,1);
  IB = 1:size(B,1);

  iter = 1;
  while true
    prev_R = R;
    prev_t = t;

    if ~isinf(max_samples)
      IA = randperm(size(A,1));
      IA = IA(1:max_samples);
      IB = randperm(size(A,1));
      IB = IB(1:max_samples);
    end

    KAB = knnsearch(KDTA,BRt(IB,:),'K',1);
    KBA = knnsearch(KDTB,bsxfun(@minus,A(IA,:),t)*R','K',1);

    [R,t] = fit_rigid([A(KAB,:);A(IA,:)],[B(IB,:);B(KBA,:)]);
    % This builtin is slower...
    %[~,~,T] = procrustes([A(KAB,:);A(IA,:)],[B(IB,:);B(KBA,:)], ...
    %  'Scaling',false,'Reflection',false);
    %R = T.T;
    %t = T.c;
    BRt = bsxfun(@plus,B*R,t);
    if iter >1 && norm(R-prev_R)<min_R_norm && norm(t-prev_t)<min_t_norm*bbd
      break;
    end
    if ~isempty(tsurf_t) && mod(iter,5)==0
      set(tsurf_t,'Vertices',BRt);
      drawnow;
      figgif('bigsig-cat-icup.gif');
      if iter == 1
      end
    end
    iter = iter + 1;
    if iter > max_iter
      if ~quiet
        warning('Max iterations (%d) exceeded without convergence',max_iter);
      end
      break;
    end
  end

  if nargout >= 4
    %KAB = knnsearch(KDTA,BRt(IB,:),'K',1);
    %KBA = knnsearch(KDTB,bsxfun(@minus,A(IA,:),t)*R','K',1);
    e = sum(sum(([A(KAB,:);A(IA,:)]-[BRt(IB,:);BRt(KBA,:)]).^2,2));
  end

end
