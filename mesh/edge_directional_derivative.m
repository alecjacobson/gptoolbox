function [DD,E,UEVperp] = edge_directional_derivative(V,F)
  % EDGE_DIRECTIONAL_DERIVATIVE
  %
  % DD = edge_directional_derivative(V,F)
  % [DD,E,UEVperp] = edge_directional_derivative(V,F)
  %
  % Inputs:
  %   V  #V by dim list of mesh vertex positions
  %   F  #F by 3 list of triangle indices into V
  % Outputs:
  %   DD  #F*3 by #V sparse matrix representing operator to comput directional
  %     derivative with respect to each edge of each face. All edges 23 then
  %     all 31 then all 12
  %   E  #F*3 by 2 list of edge vectors matching DD and UEVperp
  %   UEVperp  #F*3 by dim  Unit vectors perpendicular to edge vectors
  %
  % Example:
  %   [DD,E,UEVperp] = edge_directional_derivative(V,F);
  %   [~,C] = on_boundary(F);
  %   EBC = 0.5*(V(E(:,1),:)+V(E(:,2),:));
  %   % Treat X-coordinate as scalar field
  %   X = V(:,1);
  %   DDX = bsxfun(@times,(DD*X),UEVperp);
  %   set(tsurf(F,V),'CData',X,fphong,'EdgeColor','none');
  %   colormap(jet(20))
  %   axis equal
  %   view(2)
  %   hold on;
  %   quiver(EBC(C,1),EBC(C,2),DDX(C,1),DDX(C,2),0.2,'k','LineWidth',1.5)
  %   hold off;

  dim = size(V,2);
  
  % gradient for each triangle
  G = grad(V,F);
  
  % edges
  E = [F(:,2) F(:,3)
       F(:,3) F(:,1)
       F(:,1) F(:,2)];
  % edge vectors
  EV = V(E(:,2),:) - V(E(:,1),:);
  % unit (normalized)
  UEV = normalizerow(EV);
  
  % unit (normalized) normals
  V(:,end+1:3) = 0;
  N = normalizerow(normals(V,F));
  V = V(:,1:dim);
  
  % pad with zeros
  N  (:,end+1:3) = 0;
  UEV(:,end+1:3) = 0;
  UEVperp = cross(repmat(N,3,1),UEV,2);
  % strip padding
  UEVperp = UEVperp(:,1:dim);
  
  DD = sparse( ...
    repmat(1:size(E,1),1,dim)', ...
    reshape(bsxfun(@plus, ...
      repmat(1:size(F,1),1,3)', ...
      (0:dim-1)*size(F,1)),3*size(F,1)*dim,1), ...
    UEVperp(:), ...
    size(E,1), ...
    size(F,1)*dim) * G;
