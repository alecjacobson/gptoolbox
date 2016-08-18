function [V,F] = equilateral_tiling(x);
  % EQUILATERAL_TILING Create a tiling of the unit hexagon using equilateral
  % triangles.
  %
  % Inputs:
  %   x  number of edges on each side of the hexgon (x=1 would produce 6
  %     triangles)
  % Outputs:
  %   V  #V by 2 list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %
  [V,F] = create_regular_grid(2*x+1,2*x+1);
  V = sqrt(8)*bsxfun(@minus,V,[0.5 0.5])*sqrt(2)/2*[-1/sqrt(3) -1;1/sqrt(3) -1];
  % Outside of hexagon
  out = find(V(:,2)>V(x+1,2)+1/(2*x+1)*0.25 | V(:,2)<-V(x+1,2)-1/(2*x+1)*0.25);
  F = F(~any(ismember(F,out),2),:);
  [V,J] = remove_unreferenced(V,F);
  F = J(F);
end
