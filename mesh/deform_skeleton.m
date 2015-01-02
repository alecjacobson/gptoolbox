function [CT,BET] = deform_skeleton(C,BE,P,dQ)
  % DEFORM_SKELETON Deform a skeleton (C,BE,P) by relative rotations dQ
  %
  % Inputs:
  %   C  #C by dim list of joint positions
  %   BE  #BE by 2 list of bone edges
  %   P  #BE list of bone parents
  %   dQ  #BE by 4 list of quaternion relative rotations
  % Outputs:
  %   CT  #BE*2 by dim list of deformed bone endpoints
  %   BET  #BE by 2 list of edge indices into CT
  %

  % number of bones
  m = size(BE,1);
  [Q,T] = forward_kinematics(C,BE,P,dQ);
  % Compute a deformed version of the skeleton
  CT = C(BE,:);
  BET = [1:m;m+(1:m)]';
  for b = 1:m
    for s = 0:1
      CT(b+s*m,:) = quatrotate(Q(b,:),CT(b+s*m,:))+T(b,:);
    end
  end
end
