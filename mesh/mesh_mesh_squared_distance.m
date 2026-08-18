function [sqrd,I,C] = mesh_mesh_squared_distance( ...
    V1,F1,V2,F2,best_so_far_sqrd)
  % MESH_MESH_SQUARED_DISTANCE Compute squared distance between two meshes.
  %
  % [sqrd,I,C] = mesh_mesh_squared_distance(V1,F1,V2,F2,best_so_far_sqrd)
  %
  % Inputs:
  %   V1  #V1 by dim list of vertex positions
  %   F1  #F1 by 3 list of triangle indices into V1
  %   V2  #V2 by dim list of vertex positions
  %   F2  #F2 by 3 list of triangle indices into V2
  %   best_so_far_sqrd  optional, ignore distances that can't improve on this
  % Outputs:
  %   sqrd  squared distance between the two meshes or empty if no distance
  %     found that's less than best_so_far_sqrd
  %   I  2 by 1 list of indices of closest triangles in F1 and F2
  %   C  2 by dim list of closest points on the two triangles
  %
  % See also: aabb_aabb_squared_distance, simplex_simplex_squared_distance
  %

  if nargin < 5
    best_so_far_sqrd = inf;
    % This is significantly faster than the matlab aabb below. And for many
    % inputs will give the best answer (point-triangle case). But it doesn't
    % really speed things up to use it for `best_so_far_sqrd` because the matlab
    % aabb will usually "rediscover" the same closest pair of triangles anyway.

    % tic;
    % % First vertex set of A vs mesh of B
    % [vm1_2,i2,c2] = point_mesh_squared_distance(V1,V2,F2);
    % [d,v1] = min(vm1_2);
    % [i2,c2] = deal(i2(v1),c2(v1,:));
    % i1 = find(any(F1==v1,2));
    % c1 = V1(v1,:);
    % assert(~isempty(i1),'No face found for vertex %d',v1);
    % I = [i1;i2];
    % C = [c1;c2];

    % [vm2_1,i1,c1] = point_mesh_squared_distance(V2,V1,F1);
    % [d2,v2] = min(vm1_2);
    % if d2 < d
    %   d = d2;
    %   [i1,c1] = deal(i1(v2),c1(v2,:));
    %   i2 = find(any(F2==v2,2));
    %   c2 = V2(v2,:);
    %   assert(~isempty(i2),'No face found for vertex %d',v2);
    %   I = [i1;i2];
    %   C = [c1;c2];
    % end
    % fprintf('%30s: %f seconds\n','point_mesh_squared_distance',toc);
    % best_so_far_sqrd = d;
  end

  tic;
  [PF1B1,PF1B2] = box_each_element(V1,F1);
  [PF2B1,PF2B2] = box_each_element(V2,F2);

  [F1B1,F1B2,F1leaf] = aabb(PF1B1,PF1B2);
  [F2B1,F2B2,F2leaf] = aabb(PF2B1,PF2B2);
  fprintf('%30s: %f seconds\n','aabb',toc);

  tic;
  [aabb_sqrd,aabb_I,aabb_C] = aabb_aabb_squared_distance( ...
    F1B1,F1B2,F1leaf, ...
    F2B1,F2B2,F2leaf, ...
    @triangle_triangle_squared_distance, ...
    best_so_far_sqrd);
  fprintf('%30s: %f seconds\n','aabb_aabb_squared_distance',toc);
  if ~isempty(aabb_sqrd) && aabb_sqrd < best_so_far_sqrd
    sqrd = aabb_sqrd;
    I = aabb_I;
    C = aabb_C;
  end

  function [d,c1,c2] = triangle_triangle_squared_distance(i,j,best_so_far_sqrd)
    [d,c1,c2] = simplex_simplex_squared_distance( ...
      V1(F1(i,:),:), ...
      V2(F2(j,:),:));
    % correct return semantics even though we didn't exploit best_so_far_sqrd
    % for performance.
    if d > best_so_far_sqrd
      d = [];
      c1 = [];
      c2 = [];
    end
  end
end


  %E1 = edges(F1);
  %E2 = edges(F2);
  %[PE1B1,PE1B2] = box_each_element(V1,E1);
  %[PE2B1,PE2B2] = box_each_element(V2,E2);
  %[E1B1,E1B2,E1leaf] = aabb(PE1B1,PE1B2);
  %[E2B1,E2B2,E2leaf] = aabb(PE2B1,PE2B2);
