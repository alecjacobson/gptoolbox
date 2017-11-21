function [V,F] = delaunayize(V,F,varargin)
  % DELAUNAYIZE Make all edges Delaunay (if possible).
  %
  % [V,F] = delaunayize(V,F,varargin)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   Optional:
  %     'Tol'  followed by tolerance for `is_delaunay`
  %     'Keep'  followed by #E by 2 list of edges _not_ to flip
  %     'SplitEdges'  followed by whether to allow edge splitting.
  %     'MaxDihedralAngle'  followed by max-dihedral angle of edges to allow
  %       flip across {inf}
  % Outputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %
  % See also: is_delaunay, split_edges, flip_edges
  %

  tol= 1e-7;
  vis = false;
  keep_E = zeros(0,2);
  split_edges = false;
  max_dihedral_angle = inf;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Tol','Visualize','Keep','SplitEdges','MaxDihedralAngle'}, ...
    {'tol','vis','keep_E','split_edges','max_dihedral_angle'});
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

  while true
    if vis
      tsurf(F,V);
      drawnow;
    end
    [D,allE]= is_delaunay(V,F,'Tol',eps,'BoundaryDefault',false);
    [~,B] = on_boundary(F);
    % Find edges to be kept
    K = reshape(ismember(sort(allE,2),keep_E,'rows'),[],3);
    if isinf(max_dihedral_angle)
      C = false(size(D));
    else
      % Edges with too much curvature are also by definition "delaunay"
      % pi is flat
      [A,G] = adjacency_dihedral_angle_matrix(V,F);
      [AI,AJ,AV] = find(A);
      [GI,GJ,GV] = find(G);
      assert(isequal(AI,GI));
      assert(isequal(AJ,GJ));
      C = full(sparse(GI,GV,abs(AV-pi),size(F,1),size(F,2)))>max_dihedral_angle;
    end
    % Kept edges, boundary edges, and curved edges are by definition "delaunay"
    D = D | K;
    D = D | B;
    D = D | C;
    % list all delaunay edges
    NDE = allE(~D,:);
    % list all unique delaunay edges
    NDuE = unique(sort(NDE,2),'rows');
    NDuKB = [];
    if split_edges
      DKB = D & (K | B);
      NDKB = allE(~DKB,:);
      NDuKB = unique(sort(NDKB,2),'rows');
      if vis
        hold on;
        plot_edges(V,NDuKB,'LineWidth',4);
        hold off;
        input('');
      end
      [V,F] = split_edges(V,F,NDuKB);
    end
    if ~isempty(NDuE)
      Vvis = [];
      if vis
        Vvis = V;
      end
      Fprev = F;
      F = flip_edges(Fprev,NDuE,'AllowNonManifold',false,'V',Vvis);
      degen = sum(F(:,1)==F(:,2) | F(:,2)==F(:,3) | F(:,3)==F(:,1));
      if degen>0
        warning('flip_edges.m created degenerate %d facets. Reverting...',degen);
        Fprev = F;
      end
      if all(F(:) == Fprev(:))
        break;
      end
    elseif isempty(NDuKB)
      break;
    end
  end
end
