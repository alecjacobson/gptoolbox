function [F,C] = bfs_orient(F,V);
  % BFS_ORIENT Consistently orient faces in orientable patches using BFS
  %
  % F = bfs_orient(F,V);
  %
  % Inputs:
  %  F  #F by 3 list of faces
  % Outputs:
  %  F  #F by 3 list of faces
  %  C  #F list of component ids
  %
  % Example:
  %  RF = F;
  %  I = rand(size(F,1),1)<0.5;
  %  RF(I,:) = fliplr(RF(I,:));
  %  [FF,C] = bfs_orient(RF);
  %  FF = orient_outward(V,F,C);
  %  I = randperm(max(C))';
  %  meshplot(V,RF,'ScalarFieldF',I(C));
  %  meshplot(V,FF,'ScalarFieldF',I(C));
  %
  % See also: orient_outward, manifold_patches
  %

  [C,A] = manifold_patches(F);
  % No self matches
  A = A-diag(diag(A));
  % loop over components
  for c = 1:max(C)
    F(C==c,:) = bfs_orient_patch(F(C==c,:),A(C==c,C==c));
  end

  function FF = bfs_orient_patch(FF,AA)
    m = size(FF,1);
    seen = false(m,1);

    es = [2 3;3 1;1 2];
    Q = 1;
    while ~isempty(Q)
      % pop
      f = Q(1);
      Q = Q(2:end);
      if seen(f)
        continue;
      end
      seen(f) = true;
      % only consider manifold-edge neighbors 
      neighbors = find(AA(:,f));
      % number of neighbors
      nn = numel(neighbors);
      % loop over edges
      for ei = 1:3
        % this directed edge
        e = FF(f,es(ei,:));
        % All directed edges of neighbors
        E = [FF(neighbors,[2 3]);FF(neighbors,[3 1]);FF(neighbors,[1 2])];
        % Neighbors with matching directed edge ...
        match = neighbors(mod(find(all(bsxfun(@eq,e,E),2))-1,nn)+1); 
        % ... must be flipped
        FF(match,:) = fliplr(FF(match,:));
        % Append neighbors to queue
        Q = [Q;neighbors];
      end

      %t = tsurf(FF(seen,:),V);
      %drawnow;
      %input('');
    end
    assert(all(seen));
  end

end
