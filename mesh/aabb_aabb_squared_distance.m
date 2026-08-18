function [sqrd,I,C] = aabb_aabb_squared_distance( ...
    A1,A2,Aleaf, ...
    B1,B2,Bleaf, ...
    leaf_leaf_squared_distance, ...
    best_so_far_sqrd)
  % AABB_AABB_SQUARED_DISTANCE Compute squared distance between two AABB trees.
  %
  % [sqrd,I,C] = aabb_aabb_squared_distance( ...
  %   A1,A2,Aleaf,B1,B2,Bleaf,leaf_leaf_squared_distance,best_so_far_sqrd)
  %
  % Inputs:
  %   A1  #A by dim list of min corners of boxes
  %   A2  #A by dim list of max corners of boxes
  %   Aleaf  #A list of leaf indices (or 0 if not a leaf)
  %   B1  #B by dim list of min corners of boxes
  %   B2  #B by dim list of max corners of boxes
  %   Bleaf  #B list of leaf indices (or 0 if not a leaf)
  %   leaf_leaf_squared_distance  function that takes primitive ids a and b on
  %     the entities stored in the A and B trees respectively and returns
  %       [sqrd,ca,cb] = leaf_leaf_squared_distance(a,b,best_so_far_sqrd)
  %     where ca and cb are the closest points on primitives a and b.
  %   best_so_far_sqrd  ignore distances that can't improve on best_so_far_sqrd
  % Outputs:
  %   sqrd  squared distance between the two AABB trees or empty if no distance
  %     found that's less than best_so_far_sqrd
  %   I  2 by 1 list of indices of closest primitives in the A and B trees
  %   C  2 by dim list of closest points on the two primitives
  %
  % See also: aabb (which builds the input eytzinger AABB trees),
  % point_mesh_squared_distance, mesh_mesh_squared_distance
  %

  if nargin < 8 || isempty(best_so_far_sqrd)
    best_so_far_sqrd = inf;
  end

  assert(size(A1,2) == size(A2,2));
  assert(size(B1,2) == size(B2,2));
  assert(size(A1,2) == size(B1,2));
  assert(size(A1,1) == numel(Aleaf));
  assert(size(A2,1) == numel(Aleaf));
  assert(size(B1,1) == numel(Bleaf));
  assert(size(B2,1) == numel(Bleaf));

  sqrd = [];
  I = [];
  C = [];

  if isempty(Aleaf) || isempty(Bleaf)
    return;
  end

  best = best_so_far_sqrd;

  root_sqrd = box_box_squared_distance(1,1);
  recursive_helper(1,1,root_sqrd);

  function d = box_box_squared_distance(ai,bi)
    % Squared distance between two axis-aligned boxes. This is a lower bound
    % on the squared distance between anything contained in the two boxes.
    delta = max( ...
      A1(ai,:) - B2(bi,:), ...
      B1(bi,:) - A2(ai,:));
    delta = max(delta,0);
    d = sum(delta.^2);
  end

  function recursive_helper(ai,bi,lower_bound)

    % This pair cannot improve the current best result.
    if ~(lower_bound < best)
      return;
    end

    aleaf = Aleaf(ai);
    bleaf = Bleaf(bi);

    % Unused locations in the Eytzinger arrays.
    if aleaf < 0 || bleaf < 0
      return;
    end

    % Leaf-leaf pair: evaluate the actual primitive distance.
    if aleaf > 0 && bleaf > 0
      [d,ca,cb] = leaf_leaf_squared_distance(aleaf,bleaf,best);
      if ~isempty(d) && d < best
        best = d;
        sqrd = d;
        I = [aleaf;bleaf];
        %C = [reshape(ca,1,[]);reshape(cb,1,[])];
        C = [ca;cb];
      end
      return;
    end

    % Construct all pairs obtained by descending whichever nodes are
    % internal.
    if aleaf == 0 && bleaf == 0
      children = [ ...
        2*ai   2*bi; ...
        2*ai   2*bi+1; ...
        2*ai+1 2*bi; ...
        2*ai+1 2*bi+1];
    elseif aleaf == 0
      children = [ ...
        2*ai   bi; ...
        2*ai+1 bi];
    else
      children = [ ...
        ai 2*bi; ...
        ai 2*bi+1];
    end

    % Compute lower bounds and visit the most promising pairs first. Finding
    % a good leaf pair early tightens `best` and improves pruning.
    D = inf(size(children,1),1);
    for c = 1:size(children,1)
      cai = children(c,1);
      cbi = children(c,2);
      if ...
          cai <= numel(Aleaf) && ...
          cbi <= numel(Bleaf) && ...
          Aleaf(cai) >= 0 && ...
          Bleaf(cbi) >= 0
        D(c) = box_box_squared_distance(cai,cbi);
      end
    end

    [D,order] = sort(D);

    for k = 1:numel(order)
      if ~(D(k) < best)
        break;
      end

      c = order(k);
      recursive_helper(children(c,1),children(c,2),D(k));

      % Nothing can improve on an exact intersection.
      if best == 0
        return;
      end
    end
  end
end
