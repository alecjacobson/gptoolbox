function Q = axisangle2quat(W,TH)
  % AXISANGLE2QUAT Convert axis angle representation of rotation to quaternion
  %
  % Q = axisangle2quat(W,TH)
  % 
  % Inputs:
  %   W  list of axes, #rotations by 3
  %   TH  list of angles, #rotations by 1
  % Outputs:
  %  Q  list of rotations stored as quaternions, one for each control
  %    #rotations by 4 (1,i,j,k) 
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  %

  % See:
  % http://en.wikipedia.org/wiki/Axis-angle_representation#Unit_Quaternions
  Q = [cos(TH/2) W.*repmat(sin(TH/2),1,3)];
  Q = Q./repmat(sqrt(sum(Q.^2,2)),1,4);
end
