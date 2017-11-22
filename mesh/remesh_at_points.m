function [FF,II] = remesh_at_points(V,F,P)
  % REMESH_AT_POINTS Given a surface mesh (V,F) and a list of points on/near the
  % surface P, subdivide triangles of (V,F) so that P are now contained in the
  % vertex set.
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indices into V
  %   P  #P by 3 list of point positions
  % Outputs:
  %   FF  #FF by 3 list of face indices into [V;P]
  %   II  #FF list of indices into F revealing birth face
  %
  % Example:
  %   P = random_points_on_mesh(V,F,size(F,1));
  %   [FF,II] = remesh_at_points(V,F,P);
  %   tsurf(FF,[V;P],'CData',II)
  %

  VV = V;
  FF = F;
  II = (1:size(F,1))';
  JJ = zeros(size(V,1),1);
  L = (1:size(P,1))';
  PP = P;
  while true
    %clf;hold on;tsurf(FF,VV,'CData',II);scatter3(PP(:,1),PP(:,2),PP(:,3));hold off;pause
    % Find closest points/faces on mesh (This is a bit redundant. Really we
    % should only look at faces that _were_ affected last round, since those
    % must have had multiple points.)
    [~,IC,C] = point_mesh_squared_distance(PP,VV,FF);
    % Only keep one closest point per face
    [I,J] = unique(IC);
    % append keepers
    n = size(VV,1);
    JJ = [JJ;L(J)];
    VV = [VV;C(J,:)];
    nCI = n+(1:numel(J))';
    unaffected = setdiff(1:size(FF,1),I);
    II = [II(unaffected);repmat(II(I),3,1)];
    FF = [ ...
      FF(unaffected,:);
      FF(I,[1 2]) nCI; ...
      FF(I,[2 3]) nCI; ...
      FF(I,[3 1]) nCI];
    % Leftovers
    L = L(setdiff(1:end,J));
    PP = P(L,:);
    if isempty(PP)
      break;
    end
  end

  %% Re-order so that output is simply VV = [V;C]
  %C = zeros(size(P));
  %C(JJ(size(V,1)+1:end),:) = VV(size(V,1)+1:end,:);
  I = JJ+size(V,1);
  I(1:size(V,1)) = 1:size(V,1);
  FF = I(FF);
  %VV = [V;C];
end
