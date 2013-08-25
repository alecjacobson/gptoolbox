function [facets] = merge_planar_patches(V,F)
  % MERGE_PLANAR_PATCHES
  %
  % [facets] = merge_planar_patches(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices
  % Outputs:
  %   facets #planar patches cell of facet outline indices
  %

  A = adjacency_dihedral_angle_matrix(V,F);
  % Adjacency matrix of nearly coplanar neighbors
  AF = A>=(pi-1e-5);
  % get connected components
  C = components(AF);

  % loop over connected components
  facets = cell(max(C),1);
  for c = 1:max(C)
    facets{c} = outline(
  end

end
