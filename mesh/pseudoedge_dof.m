function [L] = pseudoedge_dof(C,PE,D,Q)
  % PSEUDOEDGE_DOF compute automatic dof (rotations, for now) for point
  % handles in a linear blend skinning deformation. Provided are only the
  % displacements (translations) at each control point.
  % 
  % [L] = pseudoedge_dof(C,PE,D,R)
  %
  % Inputs:
  %   C  #C by dim list of control point rest positions
  %   PE #PE by 2 list of indices into C for the pseudo edges 
  %   D  #C by dim list of control point displacements 
  %   Optional:
  %     Q  #C by 4 list of existing quaternions rotations at each control,
  %       [1 i j k] form
  % Outputs:
  %   L dim by dim by #C list of linear transformations
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %

  % number of point handles
  m = size(C,1);
  % number of dimensions (2,3)
  dim = size(C,2);
  assert(size(D,1) == m);
  assert(size(D,2) == dim);

  if isempty(PE)
    % trivial case, no psuedo edges all are identity
    L = repmat(eye(dim,dim+1),[1 1 m]);
    L(:,dim+1,:) = D';
    return;
  end

  if ~exist('Q','var') || isempty(Q)
    Q = repmat([1 0 0 0],m,1);
  end

  assert(size(Q,2) == 4);
  assert(size(Q,1) == m);

  % number of psuedo edges
  npe = size(PE,1);
  assert(2 == size(PE,2));
  assert(max(PE(:)) <= m);
  assert(min(PE(:)) >= 1);
  assert(all(0~=PE(:,1)-PE(:,2)));
  


  % rest vectors of each pseudo edge
  rest_PE = C(PE(:,2),:)-C(PE(:,1),:);

  % pose positions
  pose_C = C+D;
  % pose vectors of each pseudo edge
  pose_PE = pose_C(PE(:,2),:)-pose_C(PE(:,1),:);

  % compute axis and angle for each pseudo_edge
  [w,a] = axisanglebetween(rest_PE,pose_PE,[0 0 1]);
  % convert axis angles to quaternions
  Qpe = axisangle2quat(w,a)';
  % sum each component of quaternions to average
  % use sparse matrix to sum for each point handle
  I = [reshape(repmat(PE(:,1)',4,1),4*npe,1); ... 
    reshape(repmat(PE(:,2)',4,1),4*npe,1)];
  J = [repmat(1:4,1,npe*2)]';
  Qc = sparse(I,J,repmat(Qpe(:),2,1),m,4);

  % add pseudoedge rotations to existing rotations
  % for points without pseudo-edges we want to keep original rotations
  Q(PE(:),:) = quatmultiply(Q(PE(:),:),normalizerow(Qc(PE(:),:)));
  R = quat2mat(normalizerow(Q));


  % throw away extra coordinate if 2d
  if(dim == 2)
    R = R(1:dim,1:dim,:);
  end

  Rstack = reshape(permute(R,[2 1 3]),[dim dim*m])';
  T = D + stacktimes(Rstack,-C) + C;
  L = zeros([dim dim+1 m]);
  L(:,1:dim,:) = R;
  L(:,dim+1,:) = T';

end
