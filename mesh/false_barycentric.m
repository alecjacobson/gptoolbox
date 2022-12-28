function [VV,FF,I,J] = false_barycentric(V,F)
  % FALSE_BARYCENTRIC Conduct false barycentric subdivision on (V,F)
  % 
  % [VV,FF,I,J] = false_barycentric(V,F)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices into rows of V
  % Outputs:
  %   VV  #V+#F by dim list of vertex positions
  %   FF  3*#F by 3 list of triangle indices into rows of VV
  %   I  #FF list of indices into F
  %   J  #F list of indices into F of added midpoints
  %   
  %
  % See also: dual_subdivide
  %
  n = size(V,1);
  VV = [V;barycenter(V,F)];
  K = n + (1:size(F,1))';
  FF = [F(:,[2 3]) K;F(:,[3 1]) K;F(:,[1 2]) K];
  I = repmat(1:size(F,1),1,3)';
  J = F(:);
end
