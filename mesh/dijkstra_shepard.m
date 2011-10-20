function [W,dist_H_V] = dijkstra_shepard(V,F,C,p)
  % DIJKSTRA_SHEPARD
  %
  % Compute shepard weights for a list of vertices, given a list of samples and
  % optionally a denominator power value. But use "geodesic distance" as
  % defined by dijkstra's shortest path distances
  %
  % W = dijkstra_shepard(V,C)
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  #F by 3 list of face indices
  %  C  list of control vertices
  %  p  (optional) power for denominator, scalar or list same size as C {2}
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  %
  % See also:
  %   shepard
  %

  % check if either are 3D but really all z's are 0
  V_flat = size(V,2) == 3 && (sqrt(sum(V(:,3).^2)) < 1e-10);
  C_flat = size(C,2) == 3 && (sqrt(sum(C(:,3).^2)) < 1e-10);
  % is both are essentially 2D then ignore z-coords
  if((size(C,2) == 2 || C_flat) && (size(V,2) == 2 || V_flat))
    % ignore z coordinate
    V = V(:,1:2);
    C = C(:,1:2);
  end

  % number of mesh vertices
  n = size(V, 1);

  % number of control vertices
  c = size(C,1);

  assert(size(C,2) == size(V,2));
  dim = size(C,2);

  % default p value
  if(~exist('p','var'))
    p = 2;
  end

  % number of domain vertices
  n = size(V,1);

  % if p is a scalar convert it now to a list the same size as C
  if(prod(size(p)) == 1)
    p = repmat(p,c,1);
  end

  % compute distance from every vertex in the mesh to every control vertex
  D = permute(sum((repmat(V,[1,1,c]) - ...
    permute(repmat(C,[1,1,n]),[3,2,1])).^2,2),[1,3,2]);
  % use distances to determine closest mesh vertex to each control vertex
  % Cv(i) is closest vertex in V to ith control vertex in C
  [minD,Cv] = min(D);
  % if number of unique closest mesh vertices is less than total number, then
  % we have contradictory boundary conditions
  if(~all(size(unique(Cv)) == size(Cv)))
    warning('Multiple control vertices snapped to the same domain vertex');
  end

  % build adjacency matrix with edge lengths as entries
  C = adjacency_edge_cost_matrix(V,F);
  dist_H_V = repmat(minD',[1 n]) + dijkstra(C,Cv);

  % power of each control point seen by each vertex in domain
  pp = repmat(p,1,n);
  W = 1.0./((dist_H_V).^pp);

  % Handle degenrate case that a control point is on a mesh vertex
  % snap vertices close to corners
  on_sample = dist_H_V < eps;
  W(:,any(on_sample)) = 0;
  W(on_sample) = 1;

  % normalize W
  W = W./repmat(sum(W,1),c,1);

  % we've made W transpose
  W = W';

end
