function [V,F] = sphericon(r)
  % SPHERICON  Construct a mesh of a sphericon
  %
  % [V,F] = sphericon(r)
  %
  % Input:
  %   r  radius
  % Outputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices
  %

  % Half a bicone
  sam = 100;
  [CV,CF] = create_regular_grid(sam,sam,0,0);
  CF = flip_ears(CV,CF);
  % Rotate 45°
  R = [cos(pi/4) sin(pi/4);-sin(pi/4) cos(pi/4)];
  CV = bsxfun(@minus,CV*R*2/sqrt(2),[0,1]);
  % Not uniform
  %BV = [CV(:,1) real(sqrt((1-abs(CV(:,2))).^2-CV(:,1).^2)) CV(:,2)];
  % More uniform
  R = (1-abs(CV(:,2)));
  T = ((CV(:,1)+(1-abs(CV(:,2))))./((1-abs(CV(:,2)))*2)*2-1)*pi/2;
  T(isnan(T)) = 0;
  BV = [R.*sin(T) R.*cos(T) CV(:,2)];

  % Rotate and connect
  V = [BV;BV*axisangle2matrix([0 1 0],-pi/2)*axisangle2matrix([1 0 0],-pi)];
  % Rotate everything to lie flat and scale
  V = r*V*axisangle2matrix([0,1,0],pi/4);
  F = [CF;size(BV,1)+CF];

  % Glue together
  [V,F] = clean(V,F,'MinDist',1e-8,'MinArea',0,'MinAngle',0, ...
      'SelfIntersections','ignore','SmallTriangles','remove');
  % double check that it's solid
  statistics(V,F,'Fast',true)

end
