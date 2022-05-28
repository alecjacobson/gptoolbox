function [V,F,Q] = cube(x,y,z)
  % CUBE Construct a mesh of the unit cube. Sides are ordered like sides of a
  % die (one of many dice).
  %
  % [V,F] = cube(x,y,z)
  % 
  % Inputs:
  %   x  number of vertices along x-axis {2}
  %   y  number of vertices along y-ayis {x}
  %   z  number of vertices along z-azis {y}
  % Outputs:
  %   V  x*y*z by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  %   Q  #Q by 3 list of quad indices
  %
  % 
  if nargin<1
    x = 2;
  end
  if nargin<2
    y = x;
  end
  if nargin<3
    z = y;
  end

  sam = [x y;z y;x z;x z;z y;x y];
  axes = [0 1 0;0 1 0;1 0 0;1 0 0;0 1 0;0 1 0];
  angles = [0 pi/2 pi/2 -pi/2 -pi/2 pi];
  V = [];
  F = [];
  for s = 1:6
    [CV,CF] = create_regular_grid(sam(s,1),sam(s,2),0,0);
    CV(:,3) = 0;
    R = round(axisangle2matrix(axes(s,:),angles(s)));
    F = [F;size(V,1)+CF];
    V = [V;(CV-0.5)*R+0.5];
  end
  Q = [F(1:2:end-1,[1 2]) F(2:2:end,[2 3])];
  
  % Should be able to do this procedurally
  [V,~,J] = remove_duplicate_vertices(V,1e-12);
  Q = J(Q);
  F = J(F);
  % oops inside out
  F = fliplr(F);
  Q = fliplr(Q);

end
