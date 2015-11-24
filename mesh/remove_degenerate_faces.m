function [VV,FF] = remove_degenerate_faces(V,F)
  % REMOVE_DEGENERATE_FACES Remove degenerate faces from a mesh (V,F)
  % this can make combinatorially manifold meshes into non-manifold ones, but
  % should not change the "carrier" solid of solid meshes.
  % 
  % [VV,FF] = remove_degenerate_faces(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indices into V
  % Outputs:
  %   VV  #VV by 3 list of vertex positions
  %   FF  #FF by 3 list of face indices into VV
  %

  iter=1;
  max_iter = 20;
  while true
    % Snap exactly duplicate vertices
    [V,~,J] = remove_duplicate_vertices(V,0);
    F = J(F);
    % Remove combinatorially degenerate faces
    F = F(F(:,1)~=F(:,2) & F(:,2)~=F(:,3) & F(:,3)~=F(:,1),:);
    % Remove geometrically degenerate faces
    D = doublearea(V,F)==0;
    if ~any(D)
      break;
    end
    % Find degenerate faces: since we've removed duplicate vertices these
    % should only be obtuse triangles
    FD = F(D,:);
    % Flip longest edge (opposite obtuse angle)
    allE = [FD(:,[2 3]);FD(:,[3 1]);FD(:,[1 2])];
    l = edge_lengths(V,FD);
    [~,longest] = max(l,[],2);
    to_flip = allE((1:size(FD,1))'+size(FD,1)*(longest-1),:);
    to_flip = unique(sort(allE((1:size(FD,1))'+size(FD,1)*(longest-1),:),2),'rows');
    F = flip_edges(F,to_flip);
    if iter == max_iter
      error('Maximum iterations reached');
    end
    iter = iter + 1;
  end

  VV = V;
  FF = F;
end
