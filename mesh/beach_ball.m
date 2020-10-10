function [V,F,C] = beach_ball(long_div_6)
  % BEACH_BALL Generate a beach ball mesh with color
  %
  % [V,F,C] = beach_ball(long_div_6)
  %
  % Input:
  %   long_div_6  number of longitude edges divided by 6
  % Output:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into rows of V
  %   C  #F by 3 list of rgb colors
  %   
  long = long_div_6*6;
  [X,Y,Z] = sphere(long);
  [F,V] = surf2patch(X,Y,Z,'triangles');
  BC = barycenter(V,F);
  stripes = ceil(6*(atan2(BC(:,2),BC(:,1))/(2*pi)+0.5));
  top = abs(BC(:,3))>0.9772;
  CM = cbrewer('Set1',6);
  C = CM(stripes,:);
  C(top,:) = 1;
  [V,~,J] = remove_duplicate_vertices(V,eps);
  F = J(F);
  I = (F(:,2)~=F(:,3))&(F(:,3)~=F(:,1))&(F(:,1)~=F(:,2));
  F = F(I,:);
  C = C(I,:);

end
