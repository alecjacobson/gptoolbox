function [V,F] = revolve(PV,PE,n)
  % REVOLVE Create a surface of revolution around the y-axis
  %
  % [V,F] = revolve(PV,PE,n)
  %
  % Inputs:
  %   PV  #PV by 2 list of points in the xy-plane
  %   PE  #PE by 2 list of edge indices into PV
  %   n  number of angular samples 
  % Outputs:
  %   V  #V by 3 list of 3D vertex positions
  %   F  #F by 3 list of indices into rows of V
  %
  % Example:
  %   % surface of revolution from .svg file
  %   [P,C,I,F,S,SW,D] = readSVG_cubics('~/Dropbox/models/max.svg');
  %   P(:,2) = max(P(:,2))-P(:,2);
  %   P(:,1) = P(:,1)-520;
  %   [PR,CR] = flatten_splines(P,C,I,F,S,SW,D);
  %   [PR,~,J] = remove_duplicate_vertices(PR,1e-5);
  %   CR = J(CR);
  %   [PV,PE] = spline_to_poly(PR,CR,1);
  %   [V,F] = revolve(PV,PE,314);
  %   tsurf(F,V,'FaceColor',blue,falpha(0.7,0));axis equal;camlight;view(0,50);

  assert(size(PV,2) == 2,'Only 2D supported');

  % winding number Boolean
  PE2 = [PE;size(PV,1)+fliplr(PE)];
  PV2 = [PV;PV.*[-1 1]];
  [PV,PF] = cdt([PV2;0 max(PV(:,2));0 min(PV(:,2))],[PE2;size(PV,1)*2+[1 2]],'UseBoundingBox',1);
  W = winding_number(PV2,PE2,barycenter(PV,PF));
  O = outline(PF(abs(W)>0.5,:));
  PE = O(barycenter(PV(:,1),O)<0,:);
  [PV,~,~,PE] = remove_unreferenced(PV,PE);
  [V,F] = extrude(PV,PE,'Levels',n);
  V(V(:,3)==1,3) = 0;
  V = [V(:,1).*cos(V(:,3)*2*pi) V(:,2) -V(:,1).*sin(V(:,3)*2*pi)];
  [V,~,J] = remove_duplicate_vertices(V,0);F = J(F);
  F = F(doublearea(V,F)>0,:);

end
