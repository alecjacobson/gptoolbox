function [V,F,Q] = spiral_mesh(r)
% [V,F,Q] = spiral_mesh(r)
%
% Inputs:
%   r  discrete "radius" of spiral (diameter is odd = 2r+1)
% Outputs:
%   V  #V by 2 list of vertex positions
%   F  #F by 3 list of triangle indices into rows of V
%   Q  #Q by 4 list of quad indices into rows of V
%
% Example:
%   [V,~,Q] = spiral_mesh(6);
%   [V,Q] = catmull_clark(V,Q,'Iterations',3);
% 

  %I = spiral(10);
  % Generate a spiral image
  I = false(2*r+1,2*r+1);
  x = r+1;
  y = r+1;
  sx = 0:1;
  sy = 0;
  dir = 0;
  for iter = 1:2*r+mod(r+1,2)
    I(max(min(y+sy,end),1),max(min(x+sx,end),1)) = true;
    x = x+sx(end);
    y = y+sy(end);
    switch dir
    case 0
      sy = -(0:abs(sx(end))+1);
      sx = 0;
    case 1
      sx = -(0:abs(sy(end))+1);
      sy = 0;
    case 2
      sy = (0:abs(sx(end))+1);
      sx = 0;
    case 3
      sx = (0:abs(sy(end))+1);
      sy = 0;
    end
    dir = mod(dir+1,4);
  end
  [X,Y] = meshgrid(linspace(-1,1,size(I,2)+1),linspace(-1,1,size(I,1)+1));
  [Q,V] = surf2patch(X,Y,X*0);
  [V,~,~,Q] = remove_unreferenced(V(:,1:2),Q(I,:));
  F = [Q(:,[1 2 3]);Q(:,[1 3 4])];
end
