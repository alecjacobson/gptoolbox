function [V,F,C] = beach_ball(long_div_6)
% [V,F,C] = beach_ball(long_div_6)
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
