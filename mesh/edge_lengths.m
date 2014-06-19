function l = edge_lengths(V,F)
  % EDGE_LENGTHS Compute the edge lengths for a triangle mesh
  % 
  % l = edge_lengths(V,F)
  %
  % Inputs:
  %   V   #V by dim list of triangle indices
  %   F  #F by 3 list of triangles
  % Outputs:
  %   l  #F by 3 list of edge lengths corresponding to 23,31,12
  %
  assert(size(F,2) == 3)
  i1 = F(:,1); i2 = F(:,2); i3 = F(:,3);
  s12 = normrow(V(i2,:) - V(i1,:));
  s13 = normrow(V(i3,:) - V(i1,:));
  s23 = normrow(V(i3,:) - V(i2,:));
  l = [s23 s13 s12];
end
