function [C1,C2] = cubic_split(C,t)
  % CUBIC_SPLIT Split a cubic Bezier curve at parameter t, producing a left and
  % right half.
  %
  % [C1,C2] = cubic_split(C,t)
  %
  % Inputs:
  %   C  #4 by dim list of cubic Bézier control points
  %   t  parameter in [0,1] to split along
  % Outputs:
  %   C1  #4 by dim list of cubic Bézier control points of first half
  %   C2  #4 by dim list of cubic Bézier control points of second half
  %
  % See also: cubic_eval
  %
  C12 = (C(2,:)-C(1,:)).*t + C(1,:);
  C23 = (C(3,:)-C(2,:)).*t + C(2,:);
  C34 = (C(4,:)-C(3,:)).*t + C(3,:);
  C123 = (C23-C12).*t+C12;
  C234 = (C34-C23).*t+C23;
  C1234 = (C234-C123).*t+C123;
  C1 = [C(1,:);C12;C123;C1234];
  C2 = [C1234;C234;C34;C(4,:)];
end
