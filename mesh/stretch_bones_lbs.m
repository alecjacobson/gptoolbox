function [U] = stretch_bones_lbs(V,C,BE,new_C,WB,WE)
  % STRETCH_BONES_LBS Linear blend skinning with stretchable bones
  %
  % [U] = stretch_bones_lbs(V,C,BE,new_C,WB,WE)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   C  #C by dim list of control vertex rest positions
  %   BE  #B by 2 list of bone indices into C
  %   new_C #C by dim list of update postions for control vertices
  %   WB  #V by #B list of bone weights
  %   WE  #V by #B list of bone endpoint weights cooresponding to BE(:,2)
  % Output:
  %   U  #V by dim list of new vertex positions
  % 


  assert(size(V,2) == size(C,2));
  % number of dimensions
  dim = size(V,2);
  % number of domain vertices
  n = size(V,1);
  % number of bones
  m = size(BE,1);
  assert(size(BE,2) == 2);
  assert(max(BE(:)) <= size(C,1));
  assert(min(BE(:)) >= 1);
  assert(size(WB,1) == n);
  assert(size(WB,2) == m);
  assert(size(WE,1) == n);
  assert(size(WE,2) == m);
  assert(all(size(C) == size(new_C)));


  rest = C(BE(:,2),:) - C(BE(:,1),:);
  pose = new_C(BE(:,2),:) - new_C(BE(:,1),:);

  % compute rotations
  [w,a] = axisanglebetween(pose,rest,[ 0 0 1]);
  R = axisangle2matrix(w,a);
  if(dim == 2)
    R = R(1:dim,1:dim,:);
  end
  % compute stretch factor
  stretch = normrow(pose)./normrow(rest) - 1;
  % compute new position according to each bone
  % stack R as one tall dim*m by dim matrix
  RR = reshape(R,[dim dim*m])';
  % RV multiply vertices by each matrix and reshape so that RV is #V by dim by
  % #B
  RV = reshape((RR*(V'))',[n dim m]);
  % compute offsets rotation times inner translation (stretch)
  % repmat WE for each dimension and rearrange so bones are last index
  WEdim = permute(repmat(WE,[1 1 dim]),[1 3 2]);
  stretch_offset = stacktimes(RR,(repmat(stretch,[1 dim]).*rest));
  %stretch_offset = RR*((repmat(stretch,[1 dim]).*rest)');

  WE_stretch_offset = WEdim .* permute(repmat(stretch_offset,[1 1 n]),[3 2 1]);
  Roffset = stacktimes(RR,(-C(BE(:,1),:)));
  % translation part
  T = ...
    permute(repmat(new_C(BE(:,1),:) + Roffset,[1 1 n]),[3 2 1]) + ...
    WE_stretch_offset;
  % sum across bones
  WBdim = permute(repmat(WB,[1 1 dim]),[1 3 2]);
  U = sum( WBdim.*(T+RV),3);


  %%% convert point weights to end point weights
  %%WE1 = WC(:,BE(:,1));
  %%WE2 = WC(:,BE(:,2));
  %%WE = WE1 ./ (WE1 + WE2);
end
