function [U,UF,I] = plane_project(V,F)
  % Project each triangle to the plane
  %
  % [U,UF,I] = plane_project(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of mesh indices
  % Outputs:
  %   U  #F*3 by 2 list of triangle positions
  %   UF  #F by 3 list of mesh indices into U
  %   I  #V by #F such that I(i,j) = 1 implies U(j,:) corresponds to V(i,:)
  %
  % Examples:
  %   % Laplace-Beltrami is intrinsic
  %   LV = cotmatrix(V,F);
  %   [U,UF,I] = plane_project(V,F);
  %   LU = cotmatrix(U,UF);
  %   ILUI = I*LU*I';
  %   max(max(abs(LV-ILUI)))
  %

  nf = size(F,1);
  % edge lengths numbered same as opposite vertices
  l = [ ...
    sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
    sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
    sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
    ];
  l1 = l(:,3); l2 = l(:,1); l3 = l(:,2);

  % first corner at origin
  U = zeros(nf*3,2);
  % second corner along x axis
  U(nf+(1:nf),1) = l1;
  % third corner rotated onto plane
  U(nf*2+(1:nf),1) = (-l2.^2 + l3.^2 + l1.^2)./(2*l1);
  U(nf*2+(1:nf),2) = sqrt(l3.^2-U(nf*2+(1:nf),1).^2);
  UF = [1:nf;nf+(1:nf);2*nf+(1:nf)]';

  I = sparse(F(:),1:3*nf,1,size(V,1),nf*3);

end
