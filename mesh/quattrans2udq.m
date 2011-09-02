function UDQ = quattrans2udq(Q,T)
  % QUATTRANS2UDQ Convert a rotation stored as a quaternion and a translation
  % to a dual quaternion
  %
  % UDQ = quattrans2udq(Q,T)
  %
  % Inputs:
  %  Q  list of rotations stored as quaternions, one for each control
  %    #controls by 4 (1,i,j,k) 
  %  T  list of translations stored as vectors, one for each control
  %    #controls by 3
  % Output:
  %  UDQ  list of rigid transformations for each control point stored as dual
  %    quaternions
  %    2 by 4 #controls
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), Ladislav Kavan
  %

  % Adapted from: http://isg.cs.tcd.ie/kavanl/dq/dqconv.c
  % non-dual part (just copy Q)
  UDQ(1,:,:) = Q';
  % dual part
  UDQ(2,1,:) = -0.5*( T(:,1).*Q(:,2) + T(:,2).*Q(:,3) + T(:,3).*Q(:,4));
  UDQ(2,2,:) =  0.5*( T(:,1).*Q(:,1) + T(:,2).*Q(:,4) - T(:,3).*Q(:,3));
  UDQ(2,3,:) =  0.5*(-T(:,1).*Q(:,4) + T(:,2).*Q(:,1) + T(:,3).*Q(:,2));
  UDQ(2,4,:) =  0.5*( T(:,1).*Q(:,3) - T(:,2).*Q(:,2) + T(:,3).*Q(:,1));
end
