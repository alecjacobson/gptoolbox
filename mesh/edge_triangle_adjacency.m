function ET = edge_triangle_adjacency(F,E)
  % EDGE_TRIANGLE_ADJACENCY Build an edge-triangle adjacency matrix for a
  % triangle mesh.
  %
  % ET = edge_triangle_adjacency(F,E)
  %
  % Input:
  %   F  #F by 3  matrix of indices of vertices at triangle corners
  %   E  #E by 2  matrix of indices of vertices at edge endpoints
  % Output:
  %   ET #E by 2  map between an edge to its two indicent faces
  %               (-1 in column 2 if the edge is on the border)
  %
  % Example:
  %  E = edges(F);
  %  ET = edge_triangle_adjacency(F,E);

  VT = vertex_triangle_adjacency(F);

  ET = zeros(size(E,1),2);

  for i=1:size(E,1)
      CFi1 = find(VT(:,E(i,1)))';
      CFi2 = find(VT(:,E(i,2)))';

      temp = intersect(CFi1,CFi2);
      if size(temp,2) == 1
          temp = [temp -1];
      end

      ET(i,:) = temp;
  end

end
