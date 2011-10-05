function [NC,NE] = collapse_close_points(C,E,epsilon)
  % COLLAPSE_CLOSE_POINTS  collapse points closer than a given epsilon
  % updating given the edge list.
  %
  % [NC,NE] = collapse_close_points(C,E)
  % [NC,NE] = collapse_close_points(C,E,epsilon)
  % 
  % Inputs:
  %   C  #C by dim list of point positions
  %   E  #E by 2 list of edges
  %   or
  %   F  #F by 3 list of triangles
  %   epsilon  minium distance allowed after collapses are complete, default is to
  %     use fraction of maximum edge length
  % Outputs:
  %   NC  #NC by dim list of new point positions
  %   NE  #NE by 2 list of new edges
  %
  % Note: There is NO guarantee that the resulting mesh will be manifold.
  % 
  % Example:
  %   %% Collapse edges
  %   % point position list
  %   C = [1,0;0,0;-1,0;0,1;eps,eps;0,-0.1;0,-1];
  %   % edge list
  %   E = [1 2; 2 3; 4 5; 5 6;6 7];
  %   % plot original
  %   subplot(1,2,1);
  %   plot([C(E(:,1),1) C(E(:,2),1)]',[C(E(:,1),2) C(E(:,2),2)]', ... 
  %     '-','LineWidth',1);
  %   % collapse close edges
  %   [NC,NE] = collapse_close_points(C,E,0.2);
  %   % plot result
  %   subplot(1,2,2);
  %   plot([NC(NE(:,1),1) NC(NE(:,2),1)]',[NC(NE(:,1),2) NC(NE(:,2),2)]', ...
  %     '-','LineWidth',1);
  %
  %   %% Collapse faces
  %   % point position list
  %   C = [1,0;0,0;-1,0;0,1;eps,eps;0,-0.1;0,-1];
  %   % face list
  %   F = [1 4 5; 1 5 2; 1 2 6; 1 6 7;3 7 6; 3 6 2; 3 2 5; 3 5 4];
  %   % plot original
  %   subplot(1,2,1);
  %   tsurf(F,C);
  %   [NC,NF] = collapse_close_points(C,F,0.2);
  %   % plot result
  %   subplot(1,2,2);
  %   tsurf(NF,NC);
  %  
  %

  assert(size(E,2) == 2 || size(E,2) == 3);

  if ~exist('epsilon','var') || isempty(epsilon)
    if size(E,2) == 2
      EE = E;
    else
      EE = edges(E);
    end
    % maximum edge length
    maxD = max(sqrt(sum((C(EE(:,1),:) - C(EE(:,2),:)).^2,2)));
    epsilon = maxD/100;
  end

  % avoid sqrts
  sqr_eps = epsilon.^2;

  % make room for outputs
  NC = C;
  NE = E;

  % set previous count to phony value to enter while loop
  prev_count = size(NC,1)+1;

  % Continue collapse iterations until we're no longer collapsing anything
  while prev_count ~= size(NC,1)
    prev_count = size(NC,1);
    % repeat each point position, (each point as seen by each other point)
    NCC = repmat(NC,[1 1 size(NC,1)]);
    D = squeeze(sum((NCC - permute(NCC,[3 2 1])).^2,2));
    % add infs along diagonal so self-distances are excluded
    D = D + diag(diag(D)+inf);
    % collect minimum distances for each point and indices to that closest point
    [minD,q] = min(D);
    % indices of points
    p = 1:size(NC,1);
    % list of "snaps", that is we're going to snap snap(:,2) to snap(:,1), or
    % we're going to remove snap(:,2) and point all edges that use to be incident
    % on snap(:,2) to snap(:,1)
    snap = unique(sort([p(minD<sqr_eps)' q(minD<sqr_eps)'],2),'rows');
    % BUT! we might have a cascade of snaps (points snapping to other points that
    % are snapping...)
    f = 1;
    % continue to collapse snaps until we have only single snaps
    while any(f)
      % check if any snap-receivers are also snappers, and where they appear
      [f,r] = ismember(snap(:,1),snap(:,2));
      % tell these snappers to actually snap to their snap-receivers'
      % snap-receiver
      snap(f,1)=snap(r(f),1);
      % we may have created duplicate snaps, so only keep unique ones
      snap = unique(sort(snap,2),'rows');
      % cascades may be many levels high so continue this until there are no more
    end

    % First, remap edges to ignore snapped vertices
    % default is identity map
    remap = p;
    % insert snaps into remap
    remap(snap(:,2)) = snap(:,1);
    % remap edges
    NE = remap(NE);
    % Second, remove points from point list, and remap edges to new point list
    % Again, default is identity map
    valid = p;
    % A vertex is valid only if was not just snapped
    valid = valid(~ismember(valid,snap(:,2)));
    % only keep valid points
    NC = NC(valid,:);
    % now, remap edges to new point list
    % Again, default is identity map
    remap = p;
    % remap tells valid points where their new position is
    remap(valid) = 1:numel(valid);
    % remap edges
    NE = remap(NE);
    % we may have made duplicate edges, only keep unique ones
    NE = unique(sort(NE,2),'rows');
    % we may have made self-loops, remove self-loops
    NE = NE(NE(:,1)~=NE(:,2),:);
    if size(E,2) == 3
      % also check self loops accross other edges
      NE = NE(NE(:,2)~=NE(:,3),:);
      NE = NE(NE(:,3)~=NE(:,1),:);
    end
  end
end
