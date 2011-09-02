function [R] = quat2mat(q)
  % QUAT2MAT convert quaternions to 3d rotation matrices
  %
  % [R] = quat2mat(q)
  %
  % Input:
  %   q is an m by 4 list of normalized quaternions, [1 i j k]
  % Output:
  %   R is a 3 by 3 by m list of rotation matrices
  %
  % See:
  %   http://en.wikipedia.org/wiki/
  %   Conversion_between_quaternions_and_Euler_angles#Rotation_matrices
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %

  assert(size(q,2) == 4);

  R = zeros(3,3,size(q,1));
  R(1,1,:) = q(:,1).^2 + q(:,2).^2 - q(:,3).^2 - q(:,4).^2;
  R(1,2,:) = 2*(q(:,2).*q(:,3) - q(:,1).*q(:,4));
  R(1,3,:) = 2*(q(:,1).*q(:,3) + q(:,2).*q(:,4));
  R(2,1,:) = 2*(q(:,2).*q(:,3) + q(:,1).*q(:,4));
  R(2,2,:) = q(:,1).^2 - q(:,2).^2 + q(:,3).^2 - q(:,4).^2;
  R(2,3,:) = 2*(q(:,3).*q(:,4) - q(:,1).*q(:,2));
  R(3,1,:) = 2*(q(:,2).*q(:,4) - q(:,1).*q(:,3));
  R(3,2,:) = 2*(q(:,1).*q(:,2) + q(:,3).*q(:,4));
  R(3,3,:) = q(:,1).^2 - q(:,2).^2 - q(:,3).^2 + q(:,4).^2;

  %yy2 = 2.0 .* q(:,2) .* q(:,2);
  %xy2 = 2.0 .* q(:,1) .* q(:,2);
  %xz2 = 2.0 .* q(:,1) .* q(:,3);
  %yz2 = 2.0 .* q(:,2) .* q(:,3);
  %zz2 = 2.0 .* q(:,3) .* q(:,3);
  %wz2 = 2.0 .* q(:,4) .* q(:,3);
  %wy2 = 2.0 .* q(:,4) .* q(:,2);
  %wx2 = 2.0 .* q(:,4) .* q(:,1);
  %xx2 = 2.0 .* q(:,1) .* q(:,1);
  %R(1,1,:) = - yy2 - zz2 + 1;
  %R(1,2,:) = xy2 + wz2;
  %R(1,3,:) = xz2 - wy2;
  %R(2,1,:) = xy2 - wz2;
  %R(2,2,:) = - xx2 - zz2 + 1;
  %R(2,3,:) = yz2 + wx2;
  %R(2,1,:) = xz2 + wy2;
  %R(2,2,:) = yz2 - wx2;
  %R(2,3,:) = - xx2 - yy2 + 1;
end
