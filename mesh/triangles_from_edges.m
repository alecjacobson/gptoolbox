function F = triangles_from_edges(E)
  % TRIANGLES_FROM_EDGES Given a graph with undirected edges E, find all
  % 3-cliques (triangle)
  %
  % F = triangles_from_edges(E)
  % 
  % Inputs:
  %   E  #E by 2 list of undirected edges
  % Outputs:
  %   F  #F by 3 list of unoriented triangles
  %
  % See also: bfs_orient
  %
  n = max(E(:));
  E2V = sparse(repmat(1:size(E,1),2,1)',E,1,size(E,1),n);
  V2V = sparse([E(:,1) E(:,2)],[E(:,2) E(:,1)],1,n,n);
  % If a clique exists with this edge and some vertex then there will be exactly
  % two paths of length one from this edge to this vertex.
  % 
  % If there exists exactly two unique paths from an edge to a vertex then both
  % endpoints must be connected to the edge (no other way to get two paths).
  %
  % 3-clique iff #paths from edge to vertex == 2
  [I,J] = find(E2V*V2V==2);
  F = unique(sort([E(I,:) J],2),'rows');
end
