function [B1,B2,leaf] = aabb(PB1,PB2)
  % Inputs:
  %   PB1  #P by dim list of min corners of primitive boxes
  %   PB2  #P by dim list of max corners of primitive boxes
  % Outputs:
  %   B1  #B by dim list of min corners of boxes
  %   B2  #B by dim list of max corners of boxes
  %   leaf  #B list of leaf indices (or 0 if not a leaf)
  %
  % Example:
  %   [FB1,FB2] = box_each_element(V,F);
  %   [B1,B2,leaf] = aabb(FB1,FB2);
  %   [cV,~,cQ] = cube(2,2,2);
  %   [CV,CQ,~,CI] = repmesh(cV,cQ,B1,B2-B1);
  %   tsurf(CQ,CV,'Tets',0,falpha(0,1));
  %

  function hi = height_of_subtree(i,height)
    % Returns the height of the subtree rooted at node i in a complete binary tree
    % with height `height`.
    % The root is at height 0.
    % and the nodes are stored in depth-first order.

    if i == 1
      hi = height;
      return;
    end
    % if i > 1 then 

  end

  function recursive_helper(i,I)
    if isempty(I)
      return;
    end
    leaf(i) = 0;
    B1(i,:) = min(PB1(I,:),[],1);
    B2(i,:) = max(PB2(I,:),[],1);
    if numel(I) == 1
      leaf(i) = I(1);
      return;
    end
    [~,dir] = max(B2(i,:) - B1(i,:));
    split_value = median(PBC(I,dir));
    equal = I(PBC(I,dir) == split_value);
    left = [I(PBC(I,dir) < split_value); equal(1:floor(numel(equal)/2))];
    right = [I(PBC(I,dir) > split_value); equal(floor(numel(equal)/2)+1:end)];
    % Eytzinger's BFS layout
    left_i = 2*i;
    right_i = 2*i + 1;

    %% DFS layout
    %% decode height of i from i and total tree height
    %j = i;
    %height_i = height;
    %while j > 1
    %  if j > 2.^(height_i-1)
    %    j = j - 2^(height_i-1);
    %  else
    %    j = j - 1;
    %  end
    %  height_i = height_i - 1;
    %end
    %left_i = i + 1;
    %% size of left subtree:
    %left_s = 2^(height_i-1) - 1;
    %right_i = i + left_s + 1;

    % But, if you have to pass around the height of the full tree, you
    % might as well pass the height of the subtree instead and then the index
    % computation avoids the funny loop to decode the subtree height.
    recursive_helper(left_i,left);
    recursive_helper(right_i,right);
  end

  m = size(PB1,1);
  PBC = 0.5 * (PB1 + PB2);
  complete_m = 2^(ceil(log2(m))+1);
  height = ceil(log2(complete_m));
  leaf = repmat(-1,complete_m,1);
  B1 = nan(size(leaf,1), size(PB1, 2));
  B2 = nan(size(B1));
  I = (1:m)';

  recursive_helper(1,I);
end
