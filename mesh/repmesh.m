function [VV,FF] = repmesh(V,F,C)
  % REPMESH Repeat a mesh by translating it by each vector in C
  %
  % [VV,FF] = repmesh(V,F,C)
  % 
  % Inputs:
  %   V  #V by dim list of base mesh vertex positions
  %   F  #F by ss list of simplex indices into V
  %   C  #C by dim list of vectors
  % Outputs:
  %   VV  #VV*#C by dim list of mesh vertex positions
  %   FF  #F*#C by ss list of simplex indices
  %
  %
  % Example:
  %   [V,F] = subdivided_sphere(3);
  %   V = V*0.1;
  %   C = rand(10,3);
  %   [VV,FF] = repmesh(V,F,C);
  %   % Color based on original mesh vertex index
  %   tsurf(FF,VV,'CData',mod(1:size(VV,1),size(V,1))',fphong)
  %   % Color based on original mesh face index
  %   tsurf(FF,VV,'CData',mod(1:size(FF,1),size(F,1))')
  %   % Color based on which vector was used (per-vertex)
  %   tsurf(FF,VV,'CData',floor(((1:size(VV,1))-1)/size(V,1))',fphong)
  %   % Color based on which vector was used (per-face)
  %   tsurf(FF,VV,'CData',floor(((1:size(FF,1))-1)/size(F,1))')
  % 
  FF = reshape(F'+permute((0:size(C,1)-1)*size(V,1),[1 3 2]),3,[])';
  VV = reshape(V'+permute(C,[2 3 1]),3,size(V,1)*size(C,1))';
end
