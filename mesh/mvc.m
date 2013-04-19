function [lambda,W,tan_W] = mvc(V,C)
  % MVC - MEAN VALUE COORDINATES according to "Mean Value Coordinates" Floater
  % 2003 Equation 2.1
  %
  % [lambda,W] = mvc(V,C)
  %
  % Inputs:
  %  V  list of vertex positions
  %  C  list of polygon vertex positions (in counter-clockwise order)
  %
  % Outputs:
  %  lambda  weights, # vertices by # handles matrix of weights
  %  W  weights before normalization, # vertices by # handles matrix of weights
  %
  %

  % at least three control points
  assert(size(C,1)>2);
  VV = [];

  % get indices to corners and prev and next
  I = 1:size(C,1);
  prev_I = mod(I-1+size(C,1)-1,size(C,1))+1;
  next_I = mod(I,size(C,1))+1;

  % check if either are 3D but really all z's are 0
  V_flat = size(V,2) == 3 && (sqrt(sum(V(:,3).^2)) < 1e-10);
  C_flat = size(C,2) == 3 && (sqrt(sum(C(:,3).^2)) < 1e-10);
  % is both are essentially 2D then ignore z-coords
  if((size(C,2) == 2 || C_flat) && (size(V,2) == 2 || C_flat))
    % ignore z coordinate
    V = V(:,1:2);
    C = C(:,1:2);
  else
    % give dummy z coordinate to either mesh or poly
    if(size(V,2) == 2)
      V = [V zeros(size(V,1),1)];
    end
    if(size(C,2) == 2)
      C = [C zeros(size(C,1),1)];
    end

    % check that C is planar
    % average normal around poly corners
    normals = cross(C(next_I,:)-C(I,:),C(prev_I,:)-C(I,:),2);
    n = mean(normals);
    % normalize n
    n = n./repmat(norm(n),1,3);
    % take centroid as point on plane
    p = mean(C);
    dist_to_plane = dot(repmat(n,size(C,1),1),C-repmat(p,size(C,1),1),2);
    % check that poly is really coplanar
    assert(sqrt(sum(dist_to_plane.^2))<1e-10);
    dist_to_plane = dot(repmat(n,size(V,1),1),V-repmat(p,size(V,1),1),2);
    % check that poly is really coplanar
    if(sqrt(sum(dist_to_plane.^2))>1e-10)
      warning('Distance from V to plane of C is large...\n');
    end

    % change of basis
    b1 = C(2,:)-C(1,:);
    b2 = cross(n,b1);
    basis = [b1;b2;n];
    % normalize basis rows
    basis = (basis./repmat(sqrt(sum(basis.^2,2)),1,3));
    % change basis of rows vectors by right multiplying with inverse of matrix
    % with basis vectors as rows
    V = V/basis;
    C = C/basis;

    % Throw away coordinates in normal direction
    V = V(:,1:2);
    C = C(:,1:2);
  end

  % vectors from V to every C, where CmV(i,j,:) is the vector from domain
  % vertex j to handle i
  CmV = ...
    permute( ...
      permute(repmat(C,[1,1,size(V,1)]),[3,2,1]) - ...
      repmat(V,[1,1,size(C,1)]),[3,1,2]);
  % distance from V to every C, where dist_C_V(i,j) is the distance from domain
  % vertex j to handle i
  dist_C_V = sqrt(sum(CmV.^2,3));
  % distance from each corner in C to the next corner so that edge_length(i) 
  % is the distance from C(i,:) to C(i+1,:) defined cyclically
  edge_length = sqrt(sum((C - C([2:end 1],:)).^2,2));


  a_prev = ...
    atan2(CmV(prev_I,:,2),CmV(prev_I,:,1)) - atan2(CmV(I,:,2),CmV(I,:,1));
  a_next = a_prev(next_I,:);

  % mean value coordinates
  tan_W = (tan(a_prev/2.0) + tan(a_next/2.0));
  W = tan_W ./ dist_C_V;

  % handle degenerate cases
  EPSILON = 1e-10;
  % snap vertices close to corners
  on_corner = dist_C_V < EPSILON;
  W(:,sum(on_corner,1)==1) = 0;
  W(on_corner) = 1;

  % snap vertices close to segments
  % repmat segments lengths by number of domain vertices (edge length seen by
  % each vertex in domain) so that edge_length_V(i,:) = segment length on poly
  % from i to i+1
  edge_length_V = repmat(edge_length,1,size(V,1));
  % domain vertex j is on the segment from i to i+1 if the distances from vj to
  % pi and pi+1 are about 
  on_segment = ...
    abs((dist_C_V + dist_C_V(next_I,:) ./ edge_length_V) - 1) < EPSILON;
  % should only be on one segment (otherwise must be on a corner and we already
  % handled that
  on_segment = on_segment & repmat(sum(on_segment,1)==1,size(C,1),1);
  W(:,sum(on_segment,1)==1) = 0;
  W(on_segment) = dist_C_V(on_segment(prev_I,:));
  W(on_segment(prev_I,:)) = dist_C_V(on_segment);

  % normalize W
  lambda = W./repmat(sum(W,1),size(C,1),1);

  % we've made W transpose
  W = W';
  tan_W = tan_W';
  lambda = lambda';
end
