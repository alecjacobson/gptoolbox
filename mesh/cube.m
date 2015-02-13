function [V,F] = cube(x,y,z)
  % CUBE Construct a mesh of the unit cube. Sides are ordered like sides of a
  % die (one of many dice).
  %
  % [V,F] = cube(x,y,z)
  % 
  % Inputs:
  %   x  number of vertices along x-axis
  %   y  number of vertices along y-ayis
  %   z  number of vertices along z-azis
  % Outputs:
  %   V  x*y*z by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  %
  % 

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
  
  % Should be able to do this procedurally
  [V,~,J] = remove_duplicate_vertices(V,1e-12);
  F = J(F);

end
