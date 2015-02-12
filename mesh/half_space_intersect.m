function [VV,FF] = half_space_intersect(V,F,p,n)
  % HALF_SPACE_INTERSECT Intersect a closed mesh (V,F) with a half space
  % defined by a point on a plane the plane's normal.
  %
  % [VV,FF] = half_space_intersect(V,F,p,n)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices
  %   p  1 by 3 point on plane position
  %   n  1 by 3 plane normal vector
  % Outputs:
  %   VV  #VV by 3 list of new mesh vertex positions
  %   FF  #FF by 3 list of new mesh triangle indices
  %

  bbd = sqrt(sum((max(V)-min(V)).^2,2));
  % row vectors
  n = reshape(n,1,[]);
  p = reshape(p,1,[]);
  N = null(n)';
  CV = bsxfun(@plus,p,bbd*[1 1;1 -1;-1 -1;-1 1]*N);
  CF = [1 2 3;1 3 4];

  VCV = [V;CV];
  FCF = [F;size(V,1)+CF];
  [VV,FF,~,J,IM] = selfintersect(VCV,FCF);

  CF = FF(J>size(F,1),:);
  BC = barycenter(VV,CF);
  w = winding_number(V,F,BC);
  CF = CF(w>0.5,:);

  F = FF(J<=size(F,1),:);
  BC = barycenter(VV,F);
  F = F(sum(bsxfun(@times,bsxfun(@minus,BC,p),n),2)>=0,:);
  FF = [F;CF];
  [VV,~,J] = remove_duplicate_vertices(VV,1e-10);
  FF = J(FF);
  [VV,J] = remove_unreferenced(VV,FF);
  FF = J(FF);
  [VV,FF] = remesh_planar_patches(VV,FF,'MinSize',4,'Force',true);

end
