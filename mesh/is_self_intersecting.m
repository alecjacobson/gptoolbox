function flag = is_self_intersecting(V,O)
  % IS_SELF_INTERSECTING Determine if a set of edges in 2D has any intersections
  %
  % flag = is_self_intersecting(V,O)
  %
  % Inputs:
  %   V  #V by 2 list of vertex positions
  %   O  #O by 2 list of edge indices into rows of V
  % Outputs:
  %   flag  true iff some edge intersects another edge (not including shared
  %     vertices)
  %

  [V,~,~,O] = remove_unreferenced(V,O);
  % This is doing strictly more work since it's constructing the intersections
  % but it's faster and more robust than segment_segment_squared_distance
  flag = size(triangulate(V,O,'Flags','-c'),1)~=size(V,1);
end
