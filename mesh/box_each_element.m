function [B1,B2] = box_each_triangle(V,F)
  % corners of Axis-aligned boxes containing each triangle
  [B1,B2] = bounds(reshape(V(F,:)',size(V,2),[],size(F,2)),3);
  B1 = B1';
  B2 = B2';
  
end
