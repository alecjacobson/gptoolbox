function [PI,PC] = triangles_to_polygons(F)
  % [PI,PC] = triangles_to_polygons(F)
  PI = reshape(F',[],1);
  PC = cumsum([0;repmat(size(F,2),size(F,1),1)]);
end
