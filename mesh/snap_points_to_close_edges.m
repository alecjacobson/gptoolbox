function [NC,NE] = snap_points_to_close_edges(C,E,epsilon)
  % SNAP_POINTS_TO_CLOSE_EDGES snap points to edges closer than a given epsilon
  % breaking that edge into two edges and updating given the edge list.
  %
  % [NC,NE] = snap_points_to_close_edges(C,E)
  % [NC,NE] = snap_points_to_close_edges(C,E,epsilon)
  % 
  % Inputs:
  %   C  #C by dim list of point positions
  %   E  #E by 2 list of edges
  %   epsilon  minium distance allowed after collapses are complete, default is to
  %     use fraction of maximum edge length
  % Outputs:
  %   NC  #NC by dim list of new point positions
  %   NE  #NE by 2 list of new edges
  %
  % Example:
  %   %% Break edges
  %   % point position list
  %   C = [1,0;-1,0;0,1;eps,eps;0,-1;1,1;-1,-1];
  %   % edge list
  %   E = [1 2; 3 4; 4 5;6 7];
  %   % plot original
  %   subplot(1,2,1);
  %   plot([C(E(:,1),1) C(E(:,2),1)]',[C(E(:,1),2) C(E(:,2),2)]', ... 
  %     '-','LineWidth',1);
  %   % break edges at close points
  %   [NC,NE] = snap_points_to_close_edges(C,E,0.2);
  %   % plot result
  %   subplot(1,2,2);
  %   plot([NC(NE(:,1),1) NC(NE(:,2),1)]',[NC(NE(:,1),2) NC(NE(:,2),2)]', ...
  %     '-','LineWidth',1);

  if ~exist('epsilon','var') || isempty(epsilon)
    if size(E,2) == 2
      EE = E;
    else
      EE = edges(F);
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
  prev_count = size(NE,1)+1;

  % Continue edge break iterations until we're no longer breaking anything
  while prev_count ~= size(NE,1)
    prev_count = size(NE,1);
    % compute projection of each point to each line segment
    [T,sqrD] = project_to_lines(NC,NC(NE(:,1),:),NC(NE(:,2),:));
    % each vertex seen by each edge
    NCNE = repmat(NC,[1 1 size(NE,1)]);
    % edge start positions
    S = NC(NE(:,1),:);
    % edge destination positions
    D = NC(NE(:,2),:);
    % distance of each point to each edge start
    sqrDS = ...
      squeeze(sum((NCNE - permute(repmat(S,[1 1 size(NC,1)]),[3 2 1])).^2,2));
    % distance of each point to each edge dest
    sqrDD = ...
      squeeze(sum((NCNE - permute(repmat(D,[1 1 size(NC,1)]),[3 2 1])).^2,2));
    % replace distances to edges when point is closest to start or dest endpoints
    % respectively
    sqrD(T<0) = sqrDS(T<0);
    sqrD(T>1) = sqrDD(T>1);
    % inf-out self-distances
    sqrD( ...
      sub2ind(size(sqrD),[NE(:,1);NE(:,2)],[1:size(NE,1) 1:size(NE,1)]')) = inf;
    % for each edge find the closest point
    [minD,break_q] = min(sqrD);
    % mask telling whether closest point for each edge is close enough
    close = minD<sqr_eps;
    % don't let edges consider breaking at far points
    break_q(~close) = -1;
    % each edge claims a closest point, find edges which claim their closest
    % point "first"
    [~,first] = unique(break_q,'first');
    % default is to not break an edge
    break_e = false(1,size(NE,1));
    % only break first edge to claim closest point
    break_e(first) = true;
    % only break edges at close points
    break_e = break_e & close;
    % new edges are: old edges, first parts of broken edges, second parts of
    % broken edges
    NE = [ ...
      NE(~break_e,:) ; ...
      NE(break_e,1) break_q(break_e)'; ...
      break_q(break_e)' NE(break_e,2)];
  end
end
