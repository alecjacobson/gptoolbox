function F = three_cliques(A,varargin)
  % THREE_CLIQUES Compute all 3-cliques in a graph
  %
  % F = three_cliques(A)
  %
  % Inputs:
  %   A  #V by #V adjacency matrix
  %   Optional:
  %     'Triangulation' followed by whether to attempt to extract a minimal
  %     triangulation which matches the given edges by orienting patches and
  %     removing small patches which are not necessary to cover edges. Fills
  %     triangular holes. {false}
  % Outputs:
  %   F  #F by 3 list of face indices into V (unoriented)


  cleanup_triangulation = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Triangulation'}, ...
    {'cleanup_triangulation'});
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

  n = size(A,1);
  assert(size(A,2) == n);
  [EI,EJ] = find(triu(A,1));
  E = [EI,EJ];
  
  % make map from s → v, d→v
  S2V = sparse((1:size(E,1))',E(:,1),1,numel(E),n);
  D2V = sparse((1:size(E,1))',E(:,2),1,numel(E),n);
  A = A>0;
  S2N0 = S2V*A;
  D2N0 = D2V*A;
  S2N = S2N0-D2V;
  D2N = D2N0-S2V;
  E2O = (S2N & D2N);
  
  % unfortunately this will find each face three times.
  [J,I] = find(E2O');
  F = [E(I,:) J];
  F = remove_duplicate_simplices(F);

  if cleanup_triangulation
    F = bfs_orient(F);
    C = manifold_patches(F);
    counts = accumarray(C',1);
    [~,order] = sort(counts,'descend');
    IM(order) = 1:length(order);
    C = IM(C);
    
    allE = [F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];
    [E,~,EMAP] = unique(sort(allE,2),'rows');
    EMAP = reshape(EMAP,[],3);
    
    covered = false(size(E,1),1);
    keep = false(size(F,1),1);
    
    for c = 1:numel(counts)
      % Any uncovered edge in this patch
      Ic = find(C==c);
      uedges = EMAP(Ic,:);
    
      not_all_covered = ~all(covered(uedges),'all');
    
      if not_all_covered
        keep(Ic) = true;
        covered(uedges) = true;
      end
    end
    F = F(keep,:);
    
    F = bfs_orient(F);
    
    O = outline(F);
    CO = connected_components(O);
    counts = accumarray(CO',1);
    triangle_boundaries = find(counts==3);
    for t = reshape(triangle_boundaries,1,[])
      Ft = find(CO==t);
      if ismember(Ft(:,1:2),O,'rows')
        Ft = fliplr(Ft);
      end
      F = [F;Ft];
    end
  end

end
