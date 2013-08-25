function [VV,FF] = extrude(V,F)
  % EXTRUDE Extrude a 2d mesh in the z direction by 1, connecting boundaries apropriately 
  %
  % Inputs:
  %  V  #V by 2 list of 2d vertex positions
  %  F  #F by 3 list of triangle indices into V
  % Outputs:
  %  VV #VV by 3 list of 3d vertex positions
  %  FF  #FF by 3 list of triangle indices into VV
  %

  % Copy top and bottom
  VV = [V 1+0*V(:,1);V 0*V(:,1)];
  FF = [F;size(V,1)+fliplr(F)];

  % connect boundaries
  O = outline(F);
  FO = [size(V,1)+O(:,1) O(:,[2 1]); O(:,2) size(V,1)+O(:,1:2)];

  FF = [FF;FO];

end
