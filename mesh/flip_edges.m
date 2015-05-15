function [FF,I,l] = flip_edges(F,E,varargin)
  % FLIP_EDGES  Flip edges (E) of a mesh (F). Mesh should contain edges E and
  % be edge-manifold along each edges in E. Edges will be flipped according to
  % the order of rows in E: **order** matters. Flipping coincident edges
  % requires care on the callers part (e.g. to avoid creating degenerate
  % triangles).
  %
  % Inputs:
  %   F  #F by 3 list of triangle indices into some list of vertices
  %   E  #E by 2 list of undirected edges.
  %   Optional:
  %     'AllowNonManifold' followed by whether to only flip edges if the new
  %       edge does not already exist. {false}
  %     'SideLengths' followed by:
  %       l  #F by 3 list of side lengths corresponding to edges 23 31 12
  %     'MaxIter' maximum iterations of disconnected edge components
  % Outputs:
  %   FF  #F by 3 list of new facets
  %   I  #F*3 list of indices such that FF = reshape(F(I),size(F))
  %   l  #F by 3 list of side lengths corresponding to edges 23 31 12
  %

  % default values
  l = [];
  max_iter = inf;
  allow_nm = false;
  V = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'SideLengths','MaxIter','V','AllowNonManifold'}, ...
    {'l','max_iter','V','allow_nm'});
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

  if ~isempty(l)
    assert(all(size(F)==size(l)),'|Edge lengths| should match |F|');
  end

  % number of vertices
  nv = max(F(:));
  nf = size(F,1);
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
  I = 1:numel(F);

  iter = 1;
  while ~isempty(E)
    % Find one edge per connected component of edges (ideally we would find a
    % "perfect matching")
    % Should be done in reduced graph of just vertices incident on E
    ne = size(E,1);
    C = connected_components(E);
    C = C(E(:,1));
    % Use reverse order so we can take max
    E2C = sparse(1:ne,C,ne-(1:ne)+1,ne,max(C));
    mI = max(E2C,[],1);
    mI = mI(mI>0);
    % reverse
    mI = ne - mI + 1;
    EE = E(mI,:);

    if ~isempty(V)
      tsurf(FF,V,'FaceColor','r');
      hold on;
      plot_edges(V,EE,'LineWidth',3);
      hold off;
      axis equal;
      input('');
    end

    % Pop those from E
    E = E(setdiff(1:end,mI),:);

    % Need to check that resulting edge is manifold.

    % EE now has at most one edge incident on a vertex
    allE = [FF(:,[2 3]);FF(:,[3 1]);FF(:,[1 2])];
    [found,f_left] = ismember(EE,allE,'rows');
    assert(all(found),'Edges should appear in facets');
    f_left_c = floor((f_left-1)/nf)+1;
    f_left = mod(f_left-1,nf)+1;
    [found,f_right] = ismember(fliplr(EE),allE,'rows');
    if any(~found) && ~isempty(V)
      tsurf(F,V);
      hold on;
      plot_edges(V,EE(~found,:),'r','LineWidth',4);
      hold off;
    end
    assert(all(found),'Edges should appear in facets in both directions');
    f_right_c = floor((f_right-1)/nf)+1;
    f_right = mod(f_right-1,nf)+1;

    %
    %  e1----fj
    %  | \    |
    %  |  \   |
    %  |   \  |
    %  |    \ |
    % fi-----e2
    %

    % Should be safe to read and write from FF alternating, since edges should
    % not touch.
    fJ_left = sub2ind(size(F),f_right,f_right_c);
    fJ_right = sub2ind(size(F),f_left,f_left_c);

    new_EE = sort([FF(fJ_left) FF(fJ_right)],2);
    if ~allow_nm
      [found] = ismember(new_EE,uE,'rows');
      if any(found)
        warning('Ignoring edge flips that would create non-manifold edges');
        new_EE = new_EE(~found,:);
        f_left_c = f_left_c(~found);
        f_left = f_left(~found);
        f_right = f_right(~found);
        f_right_c = f_right_c(~found);
        fJ_left = fJ_left(~found);
        fJ_right = fJ_right(~found);
      end
    end
    fI_left = sub2ind(size(F),f_left,mod(f_left_c+1,3)+1);    % +2
    fI_right = sub2ind(size(F),f_right,mod(f_right_c+1,3)+1); % +2

    if nargout >= 3
      % New edge lengths via [Fisher et al. 2007]
      le = sub2ind(size(F),f_left,mod(f_left_c-1+0,3)+1);    % +0
      la = sub2ind(size(F),f_left,mod(f_left_c-1+1,3)+1);    % +1
      lb = sub2ind(size(F),f_left,mod(f_left_c-1+2,3)+1);    % +2
      re = sub2ind(size(F),f_right,mod(f_right_c-1+0,3)+1);    % +0
      rc = sub2ind(size(F),f_right,mod(f_right_c-1+1,3)+1);    % +1
      rd = sub2ind(size(F),f_right,mod(f_right_c-1+2,3)+1);    % +2
      a = l(la);
      b = l(lb);
      e = l(le); % == l(fI_right_e)
      c = l(rc);
      d = l(rd);
      % tan(alpha/2)
      tan_a_2 = sqrt(((a-b+e).*(a+b-e))./((a+b+e).*(-a+b+e)));
      tan_d_2 = sqrt(((d+e-c).*(d-e+c))./((d+e+c).*(-d+e+c)));
      % tan((alpha+delta)/2)
      tan_a_d_2 = (tan_a_2 + tan_d_2)./(1-tan_a_2.*tan_d_2);
      % cos(alpha+delta)
      cos_a_d = (1-tan_a_d_2.^2)./(1+tan_a_d_2.^2);
      f = sqrt(b.^2 + c.^2 - 2.*b.*c.*cos_a_d);

      % Need to swap, so be lazy and copy
      old_l = l;
      % vb --> ve, la --> f,  le --> lc
      l(la) = f;
      l(le) = old_l(rc);
      % vd --> ve, lc --> f, re --> la
      l(rc) = f;
      l(re) = old_l(la);
    end

    I(fI_left) = fJ_left;
    FF(fI_left) = FF(fJ_left); 
    I(fI_right) = fJ_right;
    FF(fI_right) = FF(fJ_right);

    if ~isempty(V)
      tsurf(FF,V,'FaceColor','r');
      hold on;
      plot_edges(V,new_EE,'LineWidth',3);
      hold off;
      axis equal;
      input('');
    end

    if iter >= max_iter
      break;
    end
    iter = iter + 1;
  end


end
