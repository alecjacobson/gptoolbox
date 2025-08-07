function [V,F,I] = capsule(n,a)
  % [V,F,I] = capsule(n,a)
  %
  % Inputs:
  %   n  number of faces along latitude of each hemisphere
  %   a  length of cylindrical portion
  % Outputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into rows of V
  %   I  #F list I(f) = 1,2 or 3 indicating bottom hemisphere, cylindrical part,
  %     top hemisphere, respectively.
  %
  % See also: stadium
  R = 1;
  s = 1;
  s = ceil( (2*n)/(2*pi*R)* a);
  [CV,CF] = cylinder_mesh(R,2*n,'Stacks',s);
  CV = CV.*[1 1 a];
  [X,Y,Z] = sphere(2*n);
  [SF,SV] = surf2patch(X,Y,Z,'triangles');
  [SV,~,J] = remove_duplicate_vertices(SV,1e-7);
  SF=J(SF);
  SF = SF(SF(:,1)~=SF(:,2) & SF(:,2)~=SF(:,3) & SF(:,3)~=SF(:,1),:);
  SBC = barycenter(SV,SF);
  T = SBC(:,3)>0;
  TV = SV+[0 0 a];
  TF = SF(T,:);
  BV = SV;
  BF = SF(~T,:);
  V = [BV;CV;TV];
  F = [BF;size(BV,1)+[CF;size(CV,1)+TF]];
  [V,~,J] = remove_duplicate_vertices(V,1e-7);
  F=J(F);
  [V,~,~,F] = remove_unreferenced(V,F);
  I = [repmat(1,size(BF,1),1);repmat(2,size(CF,1),1);repmat(3,size(TF,1),1)];
end
