function [U] = dualquatlbs(V,DQ,W)
  % DUALQUATLBS Compute dual quaternions linear blend skinning deformation of
  % vertices V, using rigid transformations stored as dual quaternions, DQ, at
  % some control points, propogated to the mesh using weights W.
  %
  % [U] = dualquatlbs(V,DQ,W)
  % 
  % Inputs:
  %  V  list of vertex positions
  %  DQ  list of rigid transformations for each control point stored as dual
  %    quaternions
  %    2 by 4 by #controls
  %  W  weights, # vertices by # handles matrix of weights
  % Output:
  %  U  list of new vertex positions
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: lbs
  %

  % pad V with zeros if not 3D
  if(size(V,2) ~= 3)
    was_2d = true;
    V = [V zeros(size(V,1),1)];
  else
    was_2d = false;
  end

  % number of control points
  m = size(DQ,3);
  % should be same in W
  assert(m == size(W,2));

  % number of domain vertices
  n = size(V,1);

  % See algorithm 1 in "Geometric skinning with approximate dual quaternion
  % blending" by Kavan et al

  % compute weighted combination of DQs for each domain vertex
  % DQs seen by every vertex in domain, 2 by 4 by m by n
  WDQ = permute(repmat(W',[1,1,2,4]),[3 4 1 2]) .* repmat(DQ,[1,1,1,n]);
  % sum of weighted dual quaternions, 2 by 4 by n
  VDQ = permute(sum(WDQ,3),[1 2 4 3]);
  % regular part, n by 4
  VDQ1 = permute(VDQ(1,:,:),[3 2 1]);
  % dual part, n by 4
  VDQ2 = permute(VDQ(2,:,:),[3 2 1]);
  %VDQ2 = permute(sum(WDQ(2,:,:,:),3),[4 2 1 3]);
  len = repmat(sqrt(sum(VDQ1.^2,2)),1,4);
  VDQ1 = VDQ1./len;
  VDQ2 = VDQ2./len;

  U = ...
    V + ...
    2*cross( ...
      VDQ1(:,2:4), ...
      cross(VDQ1(:,2:4),V,2)+repmat(VDQ1(:,1),1,3).*V,2) + ...
    2*( ...
      repmat(VDQ1(:,1),1,3).*VDQ2(:,2:4) - ...
      repmat(VDQ2(:,1),1,3).*VDQ1(:,2:4) + ...
      cross(VDQ1(:,2:4),VDQ2(:,2:4),2));

  %VDQ = zeros([2 size(VDQ1)]);
  %VDQ(1,:,:) = VDQ1;
  %VDQ(2,:,:) = VDQ2;
  %% convert DQs into transformations, multiply against vertices
  %% convert rotation to quaternion and translation
  %[Q,T] = udq2quattrans(VDQ);
  %% convert quaterion
  %R = quat2mat(Q);

  % unpad U if V was 2D
  if(was_2d)
    U = U(:,1:2);
  end
end
