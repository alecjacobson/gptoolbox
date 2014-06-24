function [P,PI] = farthest_points(V,k)
  % FARTHEST_POINTS Use an iterative heuristic to sample a discrete set of
  % points so that minimum pairwise distances for each point are maximized:
  %
  % maximize ∑_i min_j(‖pi-pj‖)
  %
  % Inputs:
  %   V  #V by dim list of points in euclidean space
  %   k  number of points to output
  % Outputs:
  %   P  k by dim list of farthest points sampled from V
  %   PI  k list of indices so that P = V(PI,:)
  %

  PI = ceil(rand(k,1)*size(V,1));

  scatter3(V(:,1),V(:,2),V(:,3),'.b');
  hold on;
  P = V(PI,:);
  s = scatter3(P(:,1),P(:,2),P(:,3),'or','SizeData',100,'LineWidth',5);
  hold off;
  view(2);
  axis equal;
  
  max_iter = 100;
  iter = 1;
  while true
    change = false;
    for pi = 1:numel(PI)
      old_PI_pi = PI(pi);
      % other points
      O = V(setdiff(PI,PI(pi)),:);
      % If this is slow try switch to knnsearch
      [D,I] = pdist2(O,V,'euclidean','Smallest',1);
      [~,PI(pi)] = max(D);
    end
    change = change || (old_PI_pi ~= PI(pi));
    iter = iter+1;
    if iter>max_iter
      warning('Reached max iterations (%d) without convergence',max_iter);
      break
    end
    if ~change 
      break;
    end
    % P = V(PI,:);
    % set(s, 'XData',P(:,1), 'YData',P(:,2), 'ZData',P(:,3));
    % drawnow;
  end

  P = V(PI,:);

end
