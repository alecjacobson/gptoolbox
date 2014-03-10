function [new_TN,mani,correct] = fixNEIGH(TT,TN)
  % FIXNEIGH Tetgen outputs neighbor indices inconsistently. This fixes them so
  % that new_TN(i,k) = j means that tet i shares the face opposite its kth
  % vertex with tet j.
  %
  % [new_TN] = fixNEIGH(TT,TN)
  %
  % Inputs:
  %   TT  #TT by dim+1  list of simplex indices
  %   TN  #TT by dim+1  list of tet neighbors, so that T(i,k) = j means tet i
  %     shares *some* face with tet j
  % Outputs:
  %   new_TN  #TT by dim+1 list of tet neighbors, new_TN(i,k) = j means that
  %     tet i shares the face opposite its kth vertex with tet j.
  %   mani  whether mesh is face-manifold, new_TN <- []
  %   correct  whether TN corresponds to face-neighbors, new_TN <- []
  %
  %



  mani = true;
  correct = true;
  new_TN = [];

  % Get edges
  TE = [repmat(1:size(TT,1),1,size(TT,2))' TN(:)];
  % verify edges are symmetric
  A = sparse(TE(TE(:,2)>0,1),TE(TE(:,2)>0,2),1);
  assert(size(A,1)==size(A,2));
  assert(max(max(A'-A))==0);

  %% Need to fix?
  % Interface triangles
  TTI = [TT(:,[2 3 4]); TT(:,[3 4 1]); TT(:,[4 1 2]); TT(:,[1 2 3])];
  A = sparse(TE(TE(:,2)>0,1),TE(TE(:,2)>0,2),sum(TTI(TE(:,2)>0,:),2));
  already_consistent = full(max(max(abs(A-A'))))==0;
  if already_consistent
    warning('Seems already consistent');
  end

  % Get rid of boundary edges
  TE = TE(TE(:,2)>0,:);
  TNI = zeros(size(TE,1),1);
  for c = 1:size(TT,2)
    [TNIc] = all(~bsxfun(@eq,TT(TE(:,1),c),TT(TE(:,2),:)),2);
    % Assert face-manifoldness
    if ~all(TNI(TNIc) == 0)
      mani = false;
    end
    TNI(TNIc) = c;
  end
  % Assert that edges correspond to face sharers
  if ~all(TNI>0)
    correct = false;
  end

  if ~mani || ~correct
    return;
  end

  % Rebuild TN
  new_TN = full(sparse(TE(:,1),TNI,TE(:,2),size(TN,1),size(TN,2)));
  % boundary values consistent with input
  new_TN(new_TN==0) = min(min(TN(:)),0);
end
