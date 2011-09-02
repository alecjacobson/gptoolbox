function [V,F] = cat_meshes(V1,F1,V2,F2)
  % CAT_MESHES concatenate two meshes
  %
  % [V,F] = cat_meshes(F1,V1,F2,V2)
  %
  V = [V1 ; V2];
  F = [F1 ; F2+size(V1,1)];
end

