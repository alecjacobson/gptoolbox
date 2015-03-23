function [D,S,C] = signed_distance_direction(P,V,F)
  % SIGNED_DISTANCE_DIRECTION Compute a direction which decreases signed
  % distance with respect to an input mesh (V,F) at a set of points P.
  %
  % [D,S,C] = signed_distance_direction(P,V,F)
  %
  % Inputs:
  %   P  #P by dim list of query points
  %   V  #V by dim list of mesh vertex positions
  %   F  #F by dim+1 list of mesh indices into V
  % Outputs:
  %   D  #P by dim list of normalized directions
  %   S  #P signed distances
  %   C  #P by dim list of closest points
  %

  dim = size(F,2);
  switch dim
  case 2 
    % Facets are really edges
    E = F;
    % O(n*m) too slow...
    [T,sqrD] = project_to_lines(P,V(E(:,1),:),V(E(:,2),:),'Segments',true);
    % snap to line segment
    [~,J] = min(sqrD,[],2);
    T = T(sub2ind(size(T),1:size(P,1),J'))';
    C = V(E(J,1),:) + bsxfun(@times,T,(V(E(J,2),:)-V(E(J,1),:)));
    % sign: bug in winding number...
    s = -2*(-winding_number(V,E,P))+1;
    vec = C-P;
    % signed distance direction
    D = bsxfun(@times,s,normalizerow(vec));
    % signed distance
    S = s.*normrow(D);
  case 3 
    [S,I,C,N] = signed_distance(P,V,F,'SignedDistanceType','pseudonormal');
    D = normalizerow(C-P);
    D = bsxfun(@times,sign(S),D);
    min_dist = 1e-5;
    too_close = abs(S) < min_dist;
    % Normals always point outside regardless of the eval point.
    D(too_close,:) = -N(too_close,:);
  end

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
  %S = (1-2*w);
  %D = bsxfun(@times,S,D);
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
