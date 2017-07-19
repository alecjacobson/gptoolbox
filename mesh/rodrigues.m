function [R,W] = rodrigues(N,phi)
  % RODRIGUES Construct the rotation matrix for rotating vectors around each
  % given vector.
  %
  % Inputs:
  %   N  #N by 3 list of 3D vectors
  %   phi  angle to rotate by.
  % Outputs:
  %   R  3#N by 3#N "rotation" matrix
  %   W  3#N by 3#N cross product matrix
  % 
  % Example:
  %   % Given a scalar function Z on a mesh (V,F)
  %   G = grad(V,F);
  %   N = normalizerow(normals(V,F));
  %   X = reshape(G*Z,[],3);
  %   Y = reshape(rodrigues(N,pi/2)*G*Z,[],3);
  %   BC = barycenter(V,F);
  %   tsurf(F,V,'CData',Z,fphong,'EdgeColor','none');
  %   hold on;
  %   quiver3(BC(:,1),BC(:,2),BC(:,3),X(:,1),X(:,2),X(:,3))
  %   quiver3(BC(:,1),BC(:,2),BC(:,3),Y(:,1),Y(:,2),Y(:,3))
  %   hold off;
  %   

  m = size(N,1);
  % Either uniform angle or per-vector
  if numel(phi)>1
    assert(numel(phi) == m);
    phi = diag(sparse(phi));
  end
  % https://math.stackexchange.com/a/142831
  %
  W = sparse( ...
    [1 2 0 2 0 1]*m+(1:m)', ...
    [0 0 1 1 2 2]*m+(1:m)', ...
    [N(:,3) -N(:,2) -N(:,3) N(:,1) N(:,2) -N(:,1)], ...
    3*m,3*m);
  if isempty(phi)
    R = [];
  else
    R = speye(3*m,3*m) + sin(phi)*W + (2*sin(phi/2)^2)*W*W;
  end
end
