function [E,G] = surface_graph(P,V,F)
  % SURFACE_GRAPH Connect a given set of points P lying on a surface mesh (V,F)
  % with edges E. This *approximates* constructing a generalization of a
  % geodesically weights maximal planar graph.
  %
  % E = surface_graph(P,V,F)
  % [E,G] = surface_graph(P,V,F)
  %
  % Inputs:
  %   P  #P by 3 list input points
  %   V  #V by 3 list of mesh veretex positions
  %   F  #F by 3 list of triangle indices into V
  % Outputs:
  %   E  #E by 2 list of edge indices into P
  %   G  #G by 3 list of triangle indices into P
  %  

  % remove any unreferenced vertices (common if F is just boundary_faces(T) for
  % some tet mesh (V,F)
  [V,I] = remove_unreferenced(V,F);
  F = I(F);
  % Only consider points that are actually on the surface (this is overkill if
  % the points already are, but Rahul's dataOut.Node included points deep inside
  % the object)
  P = P(point_mesh_squared_distance(P,V,F)<eps('single'),:);
  % Upsample the surface mesh. This is necessary for the PDE coming later. It
  % could be skipped if the points P are much farther apart than the resolution
  % of (V,F). You could also play with the number of iterations or try to do it
  % adaptively. These next to steps are basically a hack for "remesh the surface
  % to have nice resolution and include the sample points".
  [VV,FF] = upsample(V,F,'Iterations',2);
  % Add points into mesh
  FF = remesh_at_points(VV,FF,P);
  VV = [VV;P];
  % Still might be duplicates
  [VV,~,J] = remove_duplicate_vertices(VV,eps);
  FF = J(FF);
  %% We really ought to have a delaunay mesh, but this is slow and can probably
  %% be avoided.
  %[VV,FF] = delaunayize(VV,FF,'Keep',edges(F),'SplitEdges',false);
  % We're going to solve a minimization:
  %
  % min_W trace(W'*L*W) subject to W(P,:) = Id;
  % 
  % Essentially this is computing "harmonic coordinates", which is a cheapshot
  % "geodesic Shepard weighting" which itself is a cheapshot Geodesic Voronoi
  % diagram.
  [b,bc] = boundary_conditions(VV,FF,P);
  Wh = kharmonic(VV,FF,b,bc);
  % We will assign each vertex of the surface mesh to the point in P with the
  % largest "weight" value.
  [~,I] = max(Wh,[],2);
  % Finally, if two neighboring mesh vertices i and j have different assigned
  % points p and q then connect p and q with an edge.
  E = unique(sort(I(edges(FF)),2),'rows');
  if nargout > 1
    % Extract all 3-cycles (triangles) from this graph (nice for visualizing
    % since otherwise you'll see all of the far away edges through the gaps.
    G = triangles_from_edges(E);
  end
end
