function [D] = signed_distance_direction(P,V,F)
  % SIGNED_DISTANCE_DIRECTION Compute a direction which decreases signed
  % distance with respect to an input mesh (V,F) at a set of points P.
  %
  % D = signed_distance_direction(P,V,F)
  %
  % Inputs:
  %   P  #P by dim list of query points
  %   V  #V by dim list of mesh vertex positions
  %   F  #F by dim+1 list of mesh indices into V
  % Outputs:
  %   D  #P by dim list of normalized directions
  %

  [S,I,C,N] = signed_distance(P,V,F);
  D = normalizerow(C-P);
  min_dist = 1e-5;
  too_close = abs(S) < min_dist;
  D(too_close,:) = -N(too_close,:);

  %% Find closest points
  %[sqrD,I,C] = point_mesh_squared_distance(P,V,F);
  %% Compute barycentric coordinates on closest faces
  %B = barycentric_coordinates(C,V(F(I,1),:),V(F(I,2),:),V(F(I,3),:));
  %% Direction to closest point
  %D = normalizerow(C-P);
  %% Determine which are too close to trust
  %min_sqr_dist = 1e-10;
  %too_close = sqrD < min_sqr_dist;
  %% Determine if closest point is on vertex, edge, or face AND too clost to
  %% trust direction.
  %epsilon = 1e-15;
  %on_face = (sum(B<=epsilon,2)==0) & too_close;
  %on_edge = (sum(B<=epsilon,2)==1) & too_close;
  %on_vertex = (sum(B<=epsilon,2)==2) & too_close;
  %% Determine which vertex or edge of that face
  %[~,which_vertex] = find(B(on_vertex,:)>epsilon);
  %[~,which_edge] = find(B(on_edge,:)<=epsilon);
  %% Compute normals for faces, vertices and edges
  %% Q: Are area normals OK?
  %% H: [Bærentzen 2002/2005] suggests angle weighting.
  %N_face = normalizerow(normals(V,F));
  %N_vertex = per_vertex_normals(V,F);
  %% map to corners
  %N_vertex = N_vertex(F,:);
  %[N_edge,E,EMAP] = per_edge_normals(V,F);
  %% map to directed edges
  %N_edge = N_edge(EMAP,:);
  %% This is an expensive way to find out if inside/outside a closed mesh
  %w = winding_number(V,F,P);
  %% Flip sign for interior points
  %D = bsxfun(@times,(1-2*w),D);
  %% Use inverse normals for those that are too close
  %if any(on_vertex)
  %  D(on_vertex,:) = ...
  %    -N_vertex(sub2ind(size(F),I(on_vertex),which_vertex),:);
  %end
  %if any(on_edge)
  %  D(on_edge,:) = -N_edge(sub2ind(size(F),I(on_edge),which_edge),:);
  %end
  %if any(on_face)
  %  D(on_face,:) = -N_face(I(on_face),:);
  %end
end
