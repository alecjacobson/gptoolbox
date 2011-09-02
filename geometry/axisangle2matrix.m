function R = axisangle2matrix(w,a)
  % AXISANGLE2MATRIX  Conver axis angle rotations into corresponding rotation
  % matrices
  %
  % R = axisangle2matrix(w,a)
  %
  % Inputs:
  %   w  n by 3 list of axis vectors
  %   a  n by 1 list of angles
  % Output:
  %   R  3 by 3 by n array of rotation matrices
  %

  % For now NaNs are not allowed
  assert(~any(isnan(w(:))));
  assert(size(w,1) == size(a,1));
  n = size(w,1);
  assert(size(w,2) == 3);


  % build the rotation matrix
  s = sin(a);
  c = cos(a);
  t = 1 - c;

  w = normalizerow(w);

  x = w(:,1);
  y = w(:,2);
  z = w(:,3);
  R = zeros([3,3,n]);
  R(1,1,:) = t.*x.*x + c;
  R(2,1,:) = t.*y.*x + s.*z;
  R(3,1,:) = t.*z.*x - s.*y;

  R(1,2,:) = t.*x.*y - s.*z;
  R(2,2,:) = t.*y.*y + c;
  R(3,2,:) = t.*z.*y + s.*x;

  R(1,3,:) = t.*x.*z + s.*y;
  R(2,3,:) = t.*y.*z - s.*x;
  R(3,3,:) = t.*z.*z + c;

end
