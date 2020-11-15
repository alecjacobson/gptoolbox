function [U,G,I,J] = remove_small_components(V,F,varargin)
  % REMOVE_SMALL_COMPONENTS
  %
  % [V,F,I,J] = remove_small_components(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indices into rows of V
  % Outputs:
  %   U  #U by 3 list of vertex positions
  %   G  #G by 3 list of face indices into rows of U
  %   I  #V by 1 list of indices such that: G = I(F)
  %   J  #G by 1 list of indices into F
  %

  [~,total_vol] = centroid(V,F);
  min_vol = 0.0001*total_vol;

  [~,C] = connected_components(F);
  nc = max(C);
  vol = zeros(nc,1);
  for i = 1:nc
    Fi = F(C==i,:);
    [~,vol(i)] = centroid(V,Fi);
  end
  vol

  J = find(ismember(C,find(vol>min_vol)));
  F = F(J,:);
  [U,I] = remove_unreferenced(V,F);
  G = I(F);
end
