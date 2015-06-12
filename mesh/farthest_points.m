function [P,PI] = farthest_points(V,k,varargin)
  % FARTHEST_POINTS Use an iterative heuristic to sample a discrete set of
  % points so that minimum pairwise distances for each point are maximized:
  %
  % maximize ???_i min_j(???pi-pj???)
  %
  % [P,PI] = farthest_points(V,k)
  % [P,PI] = farthest_points(V,k,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by dim list of points in euclidean space
  %   k  number of points to output
  %   Optional:
  %     'F' followed by #F by 3 list of facet indices into V
  %     'Distance' followed by one of:
  %       {'euclidean'} Euclidean distance
  %       'biharmonic'  biharmonic distance embedding
  %       'geodesic'  fast marching geodesic distance (slow)
  % Outputs:
  %   P  k by dim list of farthest points sampled from V
  %   PI  k list of indices so that P = V(PI,:)
  %
  % See also: random_points_on_mesh
  %

  vis = false;
  distance = 'euclidean';
  F = [];
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Distance','F','Visualize'},{'distance','F','vis'});
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

  PI = ceil(rand(k,1)*size(V,1));

  if vis
    scatter3(V(:,1),V(:,2),V(:,3),'.b');
    hold on;
    P = V(PI,:);
    s = scatter3(P(:,1),P(:,2),P(:,3),'or','SizeData',100,'LineWidth',5);
    hold off;
    view(2);
    axis equal;
    drawnow;
  end

  % embedding of V
  switch distance
  case 'euclidean'
    EV = V;
  case {'biharmonic','geodesic'}
    assert(~isempty(F));
    switch distance
    case 'biharmonic'
      EV = biharmonic_embedding(V,F,10);
    case 'geodesic'
      EV = V;
    end
  end
  
  [I,D] = knnsearch(EV(PI,:),EV,'K',1);
  max_iter = 100;
  iter = 1;
  while true
    change = false;
    %progressbar(iter-1,max_iter-1,30);
    for pi = 1:numel(PI)
      old_PI_pi = PI(pi);
      % other points
      others = PI([1:pi-1 pi+1:end]);
      % There should be a way to only compute new distances in places that
      % need to be updated...
      switch distance
      case {'euclidean','biharmonic'}
        O = EV(others,:);
        if isempty(I)
          Ipi = true(size(V,1),1);
          J = [1:pi-1 pi+1:k];
        else
          Ipi = I==pi;
          J = setdiff(1:numel(PI),pi);
        end
        %[D,I] = pdist2(O,EV,'euclidean','Smallest',1);
        % Much faster for large O/EV
        %[I,D] = knnsearch(O,EV,'K',1);
        [IIpi,D(Ipi)] = knnsearch(EV(PI(J),:),EV(Ipi,:),'K',1);
        I(Ipi) = J(IIpi);
        fIpi = find(Ipi);
        [~,d] = max(D(fIpi));
        PI(pi) = fIpi(d);
        Dpi = sqrt(sum(bsxfun(@minus,EV(PI(pi),:),EV).^2,2));
        Cpi = Dpi<D;
        D(Cpi) = Dpi(Cpi);
        I(Cpi) = pi;
      case 'geodesic'
        [D,~,I] = perform_fast_marching_mesh(V,F,others,struct('nb_iter_max',inf));
        [~,PI(pi)] = max(D);
      end
      change = change || (old_PI_pi ~= PI(pi));
    end
    if vis
      P = V(PI,:);
      %set(s, 'XData',P(:,1), 'YData',P(:,2), 'ZData',P(:,3));
      tsurf(F,V,'CData',I,'EdgeColor','none');
      hold on;
      scatter3(V(PI,1),V(PI,2),V(PI,3),'or','SizeData',100,'LineWidth',2);
      hold off;
      drawnow;
    end
    iter = iter+1;
    if iter>max_iter
      warning('Reached max iterations (%d) without convergence',max_iter);
      break
    end
    if ~change 
      break;
    end
  end

  P = V(PI,:);

end
