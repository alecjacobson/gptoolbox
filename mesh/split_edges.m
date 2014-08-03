function [U,FF] = split_edges(V,F,E)
  % SPLIT_EDGES Split given edges E in their barycenters within a mesh (V,F)
  %
  % Inputs:
  %   V  #V by dim list of mesh vertices
  %   F  #F by 3 list of mesh corners
  %   E  #E by 2 list of edges
  % Outputs:
  %   U  #V+#E by dim list of new mesh vertices
  %   FF #FF by 3 list of new mesh face corners into U
  %
  % Known bugs: this has not been tested for non-manifold edges in E, though it
  % should probably work.
  % 

  max_iter = inf;
  vis = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MaxIter','V'}, ...
    {'max_iter','V'});
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

  % number of vertices
  nv = max(F(:));
  assert(max(E(:))<=nv,'Edges exceed vertices in faces');

  % Direct adjacency matrix
  allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
  sortE = sort(allE,2);
  uE = unique(sortE,'rows');

  DA = sparse(sortE(:,1),sortE(:,2),1-2*(allE(:,2)>allE(:,1)),nv,nv);
  % Edge partiy sum adjacency matrix
  EP = DA - DA';
  % Symmetrize
  EP = EP + EP';
  % Non-manifold edges in query
  NM = EP(sub2ind(size(EP),E(:,1),E(:,2)))~=0;
  if any(NM)
    warning('Non-manifold, boundary edge flips ignored.');
    E = E(~NM,:);
  end

  FF = F;

  iter = 1;
  while ~isempty(E)
    % Find one edge per connected component of edges (ideally we would find a
    % "perfect matching")
    % Should be done in reduced graph of just vertices incident on E
    ne = size(E,1);
    C = components(adjacency_matrix(E));
    C = C(E(:,1));
    % Use reverse order so we can take max
    E2C = sparse(1:ne,C,ne-(1:ne)+1,ne,max(C));
    mI = max(E2C,[],1);
    mI = mI(mI>0);
    % reverse
    mI = ne - mI + 1;
    EE = E(mI,:);

    if vis
      tsurf(FF,V,'FaceColor','r');
      hold on;
      plot_edges(V,EE,'LineWidth',3);
      hold off;
      axis equal;
      input('');
    end

    % Pop those from E
    E = E(setdiff(1:end,mI),:);

    nf = size(FF,1);
    % EE now has at most one edge incident on a vertex
    allE = [FF(:,[2 3]);FF(:,[3 1]);FF(:,[1 2])];
    [found,f_left] = ismember(EE,allE,'rows');
    assert(all(found),'Edges should appear in facets');
    f_left_c = floor((f_left-1)/nf)+1;
    f_left = mod(f_left-1,nf)+1;
    [found,f_right] = ismember(fliplr(EE),allE,'rows');
    assert(all(found),'Edges should appear in facets in both directions');
    f_right_c = floor((f_right-1)/nf)+1;
    f_right = mod(f_right-1,nf)+1;

    %    c
    %   /|\
    %  / | \
    % a--d--b
    %

    BC = barycenter(V,EE);
    nv = size(V,1);
    V = [V;BC];
    old_nf = size(FF,1);
    FF = [FF; ...
      FF(f_right,:); ...
      FF(f_left,:)];
    nf = size(FF,1);
    FF(sub2ind([nf 3],f_right,mod(f_right_c+1,3)+1)) = nv+(1:size(EE,1));
    FF(sub2ind([nf 3],f_left,mod(f_left_c+1,3)+1)) = nv+(1:size(EE,1));
    FF(sub2ind([nf 3],old_nf+(1:size(EE,1))',mod(f_right_c,3)+1)) = nv+(1:size(EE,1));
    FF(sub2ind([nf 3],old_nf+size(EE,1)+(1:size(EE,1))',mod(f_left_c,3)+1)) = nv+(1:size(EE,1));

    if vis
      tsurf(FF,V,'FaceColor','r','FaceAlpha',0.5);
      axis equal;
      input('');
    end

    if iter >= max_iter
      break;
    end
    iter = iter + 1;

  end

  U = V;
end
