function [c,r,D] = sphere_tree(P,R)
  % [c,r,D] = sphere_tree(P,R)
  %
  % Computes a tree of spheres that cover the points in P with radii R.
  %
  % Inputs:
  %   P  #P by dim list of points
  %   R  #P by 1 list of radii (optional, default is zeros)
  % Outputs:
  %   c  #c by dim list of centers of spheres
  %   r  #c by 1 list of radii of spheres
  %   D  #c by 2 list of indices, where D(i,:) are the left and right child
  %      of the i-th sphere in the tree. If D(i,1) is NaN, then the i-th sphere
  %      is a leaf and D(i,2) is the index of the point that it covers.
  %



  if nargin<2
    R = zeros(size(P,1),1);
  end
  if size(P,1) == 1
    c = P;
    r = R;
    D = [nan 1];
    return;
  end

  % Compute the root tree;
  [c_root,r_root] = minimal_bounding_sphere(P,R);

  assert(all( r_root - (sqrt(sum((P-c_root).^2,2)) + R)  > -1e-6 ));

  % split into two groups
  bb = max(P) - min(P);
  % split along the largest dimension
  [bb_max,bb_dim] = max(bb);
  % sort along bb_dim
  [~,I] = sort(P(:,bb_dim));
  mid = floor(size(P,1)/2);
  I1 = I(1:mid,:);
  I2 = I(mid+1:end,:);
  % recurse
  [c1,r1,D1] = sphere_tree(P(I1,:),R(I1));
  [c2,r2,D2] = sphere_tree(P(I2,:),R(I2));
  % map indices to original indices
  D1(isnan(D1(:,1)),2) = I1(D1(isnan(D1(:,1)),2));
  D2(isnan(D2(:,1)),2) = I2(D2(isnan(D2(:,1)),2));
  D1(~isnan(D1(:,1)),:) = D1(~isnan(D1(:,1)),:) + 1;
  D2(~isnan(D2(:,1)),:) = D2(~isnan(D2(:,1)),:) + 1 + size(c1,1);

  % combine
  c = [c_root;c1;c2];
  r = [r_root;r1;r2];
  D = [1+1 1+size(c1,1)+1; D1; D2];

end
