function [I,minD,VI] = snap_points(C,V)
  % SNAP_POINTS snap list of points C to closest of another list of points V
  %
  % [I,minD,VI] = snap_points(C,V)
  % 
  % Inputs:
  %   C  #C by dim list of query point positions
  %   V  #V by dim list of data point positions
  % Outputs:
  %   I  #C list of indices into V of closest points to C
  %   minD  #C list of squared distances to closest points
  %   VI  #C by dim list of new point positions, VI = V(I,:)
  %

  assert(size(V,2) == size(C,2));

  % number of mesh vertices
  n = size(V, 1);

  % number of control vertices
  c = size(C,1);

  % compute distance from every vertex in the mesh to every control vertex
  D = permute(sum((repmat(V,[1,1,c]) - ...
    permute(repmat(C,[1,1,n]),[3,2,1])).^2,2),[1,3,2]);
  % use distances to determine closest mesh vertex to each control vertex
  % Cv(i) is closest vertex in V to ith control vertex in C
  [minD,I] = min(D);
  if(nargout == 3)
    VI = V(I,:);
  end
end

