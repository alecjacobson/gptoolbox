function [VV,FF,II,JJ] = remesh_at_points(V,F,P)
  VV = V;
  FF = F;
  II = (1:size(F,1))';
  JJ = zeros(size(V,1),1);
  L = (1:size(P,1))';
  PP = P;
  while true
    %clf;hold on;tsurf(FF,VV,'CData',II);scatter3(PP(:,1),PP(:,2),PP(:,3));hold off;pause
    % Find closest points/faces on mesh
    [~,I,C] = point_mesh_squared_distance(PP,VV,FF);
    % Only keep one closest point per face
    [I,J] = unique(I);
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
  % Now that we've meshed everything using closest points. Snap back to actually
  % input points.
  %VV(JJ>0,:) = P(JJ(JJ>0),:);
end
