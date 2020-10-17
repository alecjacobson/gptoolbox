function C = rotate_about(A,B,th)
  % Rotate vectors A about axis defined by vector B by amount th
  %
  % C = rotate_about(A,B,th)
  %
  % Inputs:
  %   A  #A by 3 list of points in space
  %   B  #A|1 by 3 list of axis vectors (not necessarily normalized)
  %   th  #A|1 list of rotation angles
  % Outputs
  %   C  #A by 3 list of positions
  %
  function r = cross2(a,b)
    % Optimizes r = cross(a,b,2), that is it computes cross products per row
    % Faster than cross if I know that I'm calling it correctly
    r =[a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
        a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
        a(:,1).*b(:,2)-a(:,2).*b(:,1)];
  end
  % http://math.stackexchange.com/a/1432182/81266
  Apara = sum(A.*B,2)./sum(B.*B,2).*B;
  Aperp = A-Apara;
  W = cross2(B,Aperp);
  X1 = cos(th)./normrow(Aperp);
  X2 = sin(th)./normrow(W);
  C = normrow(Aperp).*(X1.*Aperp + X2.*W) +Apara;
end
