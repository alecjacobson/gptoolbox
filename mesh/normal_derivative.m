function [DD,E] = normal_derivative(V,F)
  % NORMAL_DERIVATIVE Computes the directional derivative **normal** to **all**
  % (half-)edges of a triangle mesh (not just boundary edges). These are
  % integrated along the edge: they're the per-face constant gradient dot the
  % rotated edge vector (unit rotated edge vector for direction then magnitude
  % for integration).
  %
  % DD = normal_derivative(V,F)
  % [DD,E,N] = normal_derivative(V,F)
  %
  % Inputs:
  %   V  #V by dim list of mesh vertex positions
  %   F  #F by 3 list of triangle indices into V
  % Outputs:
  %   DD  #F*3 by #V sparse matrix representing operator to compute directional
  %     derivative with respect to each edge of each face. All edges 23 then all
  %     31 then all 12
  %   E  #F*3 by 2 list of edge vectors matching DD and N
  %   %N  #F*3 by dim  vectors perpendicular to edge vectors
  %

%  % Example:
%   %   [DD,E,N] = normal_derivative(V,F);
%   %   % normalize and Divide by lengths again to account for integral averaging
%   %   UN = bsxfun(@rdivide,N,normrow(N).^2);
%   %   [~,C] = on_boundary(F);
%   %   EBC = 0.5*(V(E(:,1),:)+V(E(:,2),:));
%   %   % Treat X-coordinate as scalar field
%   %   X = V(:,1);
%   %   DDX = bsxfun(@times,(DD*X),UN);
%   %   tsurf(F,V,'CData',X,fphong,'EdgeColor','none');
%   %   colormap(jet(20))
%   %   axis equal
%   %   view(2)
%   %   hold on;
%   %   quiver(EBC(C,1),EBC(C,2),DDX(C,1),DDX(C,2),0.2, ...
%   %     'k','LineWidth',2)
%   %   hold off;
%   %   title('Integral average normal derivative','FontSize',20);


  dim = size(V,2);

  ss = size(F,2);
  switch ss
  case 4
    % Elements are really tets
    T = F;
    % 23,31,12,41,42,43
    C = cotangent(V,T);
    % -1/1 point to/from vertex on "this edge"
    F = [T(:,[2 3 4]);T(:,[1 4 3]);T(:,[1 2 4]);T(:,[1 3 2])];
    m = size(T,1);
    DD = 2*sparse( ...
      reshape(bsxfun(@plus,repmat(1:m,1,6)',m*(0:4-1)),m,6*4), ...
      T(:,[2 1 3 1 4 1,3 2 4 2 1 2,4 3 1 3 2 3,1 4 2 4 3 4]), ...
      C(:,[3 3 2 2 4 4,1 1 5 5 3 3,6 6 2 2 1 1,4 4 5 5 6 6])* ...
        diag(repmat([1,-1],1,3*4)), ...
        m*4,size(V,1));
    E = F;
    %% REFERENCE VERSION:
    %% gradient for each tetrahedron
    %G = grad(V,F);
    %% edges
    %E = [F(:,[2 3 4]);F(:,[1 4 3]);F(:,[1 2 4]);F(:,[1 3 2])];
    %% un-normalized normals
    %N = normals(V,E);
    %V = V(:,1:dim);
    %DD = sparse( ...
    %  repmat(1:size(E,1),1,dim)', ...
    %  reshape(bsxfun(@plus, ...
    %    repmat(1:size(F,1),1,4)', ...
    %    (0:dim-1)*size(F,1)),4*size(F,1)*dim,1), ...
    %  N(:), ...
    %  size(E,1), ...
    %  size(F,1)*dim)/3 * G;
  case 3
    % edges
    E = [F(:,2) F(:,3)
         F(:,3) F(:,1)
         F(:,1) F(:,2)];
    % Could replace with:
    C = cotangent(V,F);
    % -1/1 point to/from vertex on "this edge"
    m = size(F,1);
    DD = 2*sparse( ...
      reshape(bsxfun(@plus,repmat(1:m,1,4)',m*(0:3-1)),m,4*3), ...
      F(:,[3 1 2 1,1 2 3 2,2 3 1 3]), ...
      C(:,[2 2 3 3,3 3 1 1,1 1 2 2])*diag(repmat([1,-1],1,2*3)), ...
      m*3,size(V,1));
    %% REFERENCE VERSION:
    %% gradient for each triangle
    %G = grad(V,F);
    %% edge vectors
    %EV = V(E(:,2),:) - V(E(:,1),:);
    %% unit (normalized) normals
    %V(:,end+1:3) = 0;
    %FN = normalizerow(normals(V,F));
    %V = V(:,1:dim);
    %% pad with zeros
    %FN  (:,end+1:3) = 0;
    %EV(:,end+1:3) = 0;
    %N = cross(EV,repmat(FN,3,1),2);
    %% strip padding
    %N = N(:,1:dim);
    %DD = sparse( ...
    %  repmat(1:size(E,1),1,dim)', ...
    %  reshape(bsxfun(@plus, ...
    %    repmat(1:size(F,1),1,3)', ...
    %    (0:dim-1)*size(F,1)),3*size(F,1)*dim,1), ...
    %  N(:), ...
    %  size(E,1), ...
    %  size(F,1)*dim) * G;
 end
end
