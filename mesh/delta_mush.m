function [U,Dl,I] = delta_mush(V,F,U0,varargin)
  % DELTA_MUSH Deform a mesh according to "Delta Mush: Smoothing Deformations
  % While Preserving Detail" [Mancewicz et al. 2014].
  %
  % U = delta_mush(V,F,U0)
  % [U,Dl,I] = delta_mush(V,F,U0,'ParameterName',parameter_value,...)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  %   U0  #V by 3 list of vertex positions of coarsely deformed mesh
  %   Optional:
  %     'Delta' followed by #V by 3 list of local delta values from smoothed
  %       (V,F) mesh (only really makes sense if also passing 'Neighbors'.
  %     'Neighbors' followed by precomputed I (see output below)
  %     'Smooth'  Followed by a function handle to a smoothing routine, e.g.
  %        @(U)conformalized_mean_curvature_flow(V,F,'MaxIter',10), but by
  %        %default: {@(U)laplacian_smooth(V,F,'cotan',[],0.1,'implicit',V,100)}
  %        default uses 10 iterations of inverse distance Laplacian smoothing
  %        as described in [Mancewicz et al. 2014].
  % Outputs:
  %   U  #V by 3 list of final deformed vertex positions
  %   I  #V list of indices of neighbor used to determine normal plane.
  %

  function US = delta_mush_smooth(V,F,U,max_iter)
    US = U;
    implicit = false;
    for iter = 1:max_iter
      if iter ==1 || implicit
        A = adjacency_edge_cost_matrix(V,F);
        [AI,AJ,AV] = find(A);
        A = sparse(AI,AJ,AV.^-1);
        A = diag(sparse(sum(A,2).^-1))*A;
      end
      US = A*US;
      if implicit
        V = A*V;
      end
    end
  end

  % default values
  I = [];
  Dl = [];
  smooth_fun = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Neighbors','Smooth','Delta'}, ...
    {'I','smooth_fun','Dl'});
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
  if isempty(smooth_fun)
    %lambda = 0.1;
    %max_iter = 100;
    %smooth_fun = ...
    %  @(U) laplacian_smooth(V,F,'cotan',[],lambda,'implicit',U,max_iter);
    smooth_fun = @(U) delta_mush_smooth(V,F,U,10);
  end

  in_frame = @(D,T,N,B) ...
    [ ...
      sum(D.*[T(:,1) N(:,1) B(:,1)],2) ...
      sum(D.*[T(:,2) N(:,2) B(:,2)],2) ...
      sum(D.*[T(:,3) N(:,3) B(:,3)],2) ...
    ];

  if isempty(Dl)
    VS = smooth_fun(V);
    [T,N,B,I] = per_vertex_frames(VS,F,'Neighbors',I);
    % Get "inverse frame"
    TNB = cat(3,T,N,B);
    ITNB = permute(TNB,[1 3 2]);
    IT = ITNB(:,:,1);
    IN = ITNB(:,:,2);
    IB = ITNB(:,:,3);

    Dg = (V-VS);
    Dl = in_frame(Dg,IT,IN,IB);
  end
    
  US = smooth_fun(U0);
  [T,N,B] = per_vertex_frames(US,F,'Neighbors',I);
  U = US+in_frame(Dl,T,N,B);

end
