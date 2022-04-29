function [c,S,E,flag] = convexness(V,F,n)
  % [c,S,E,flag] = convexness(V,F,n) Give a score of "how convex" a shape (V,F)
  % is on a scale of 0 to 1. "Weak Convex Decomposition by Lines-of-sight"
  % [Asafi et al. 2013]
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh face indices into rows of V
  %   n  number of samples to use {1000}
  % Outputs:
  %   c  convexity measure ∈[0,1]
  %   S  #S by n sample locations
  %   E  #S-choose-2 by 2 edge indices into rows of S
  %   flag #E list of flags whether line of sight is clear
  % 
  N = normalizerow(normals(V,F));
  [S,I] = random_points_on_mesh(V,F,n);
  % bump points a bit inside
  bbd = normrow(max(V)-min(V));
  S = S-N(I,:)*1e-4*bbd;
  % Reject any that are now outside
  W = winding_number(V,F,S,'Fast',true);
  S = S(abs(W)>0.5,:);
  % shoot rays for all pairs
  E = nchoosek(1:size(S,1),2);
  dir = S(E(:,2),:)-S(E(:,1),:);
  tic;
  [~,T] = ray_mesh_intersect(S(E(:,1),:),dir,V,F);
  toc
  % Find legit hits
  epsilon = 1e-8*bbd;
  bad = T>epsilon & T<1-epsilon;
  flag = ~bad;
  c = mean(flag);
end
