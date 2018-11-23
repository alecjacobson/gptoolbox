function [R,C,B] = circumradius(V,T)
  % CIRCUMRADIUS Return the circumradius of each triangle/tet element
  % 
  % R = circumradius(V,T)
  % [R,C,B] = circumradius(V,T)
  %
  % Input:
  %   V  #V by dim list of vertex positions
  %   T  #T by simplex-size list of tet indices
  % Output:
  %   R  #T by 1 list of simplex circumradii
  %   C  #T by dim list of simplex circumcenters
  %   B  #T by simplex-size list of barycentric coordinates so that: 
  %     C(t,:) = B(t,:) * V(T(t,:),:)
  %
  % Known issues:
  %   B output only supported for triangles
  %
  
  switch size(T,2)
  case 4
    d_14 = sqrt(sum((V(T(:,1),:)-V(T(:,4),:)).^2,2));
    d_23 = sqrt(sum((V(T(:,2),:)-V(T(:,3),:)).^2,2));
    d_24 = sqrt(sum((V(T(:,2),:)-V(T(:,4),:)).^2,2));
    d_34 = sqrt(sum((V(T(:,3),:)-V(T(:,4),:)).^2,2));
    d_12 = sqrt(sum((V(T(:,1),:)-V(T(:,2),:)).^2,2));
    d_13 = sqrt(sum((V(T(:,1),:)-V(T(:,3),:)).^2,2));
    R = sqrt((d_14.*d_23+d_13.*d_24-d_12.*d_34).*(d_14.*d_23-d_13.*d_24+d_12.*d_34).*(-d_14.*d_23+d_13.*d_24+d_12.*d_34).*(d_14.*d_23+d_13.*d_24+d_12.*d_34))./(sqrt(2)*sqrt(-2.*d_34.^2.*d_12.^4-2.*d_34.^4.*d_12.^2-2.*d_13.^2.*d_23.^2.*d_12.^2+2.*d_14.^2.*d_23.^2.*d_12.^2+2.*d_13.^2.*d_24.^2.*d_12.^2-2.*d_14.^2.*d_24.^2.*d_12.^2+2.*d_13.^2.*d_34.^2.*d_12.^2+2.*d_14.^2.*d_34.^2.*d_12.^2+2.*d_23.^2.*d_34.^2.*d_12.^2+2.*d_24.^2.*d_34.^2.*d_12.^2-2.*d_14.^2.*d_23.^4-2.*d_13.^2.*d_24.^4-2.*d_14.^4.*d_23.^2+2.*d_13.^2.*d_14.^2.*d_23.^2-2.*d_13.^4.*d_24.^2+2.*d_13.^2.*d_14.^2.*d_24.^2+2.*d_13.^2.*d_23.^2.*d_24.^2+2.*d_14.^2.*d_23.^2.*d_24.^2-2.*d_13.^2.*d_14.^2.*d_34.^2+2.*d_14.^2.*d_23.^2.*d_34.^2+2.*d_13.^2.*d_24.^2.*d_34.^2-2.*d_23.^2.*d_24.^2.*d_34.^2));
    assert(nargout == 1);
  case 3
    % http://www.mathworks.com/matlabcentral/fileexchange/17300-circumcircle-of-a-triangle/content/circumcircle.m

    %compute the length of sides (AB, BC and CA) of the triangle
    l = [ ...
      sqrt(sum((V(T(:,2),:)-V(T(:,3),:)).^2,2)) ...
      sqrt(sum((V(T(:,3),:)-V(T(:,1),:)).^2,2)) ...
      sqrt(sum((V(T(:,1),:)-V(T(:,2),:)).^2,2)) ...
      ];
    dblA = doublearea(V,T);
    %use formula: R=abc/(4*area) to compute the circum radius
    R = prod(l,2)./(2*dblA);
    %compute the barycentric coordinates of the circum center
    B= [ ...
      l(:,1).^2.*(-l(:,1).^2+l(:,2).^2+l(:,3).^2) ...
      l(:,2).^2.*( l(:,1).^2-l(:,2).^2+l(:,3).^2) ...
      l(:,3).^2.*( l(:,1).^2+l(:,2).^2-l(:,3).^2)];
    % normalize
    B = bsxfun(@rdivide,B,sum(B,2));
    %convert to the real coordinates
    C = zeros(size(T,1),size(V,2));
    for d = 1:size(V,2);
      C(:,d) = sum(B.*[V(T(:,1),d) V(T(:,2),d) V(T(:,3),d)],2);
    end

end
