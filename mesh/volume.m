function [v,sv] = volume(V,T)
  % VOLUME Compute volumes of tets T defined over vertices V
  %
  % v = volume(V,T)
  % 
  % Inputs:
  %   V  #V by dim>=3 list of vertex positions
  %   T  #T by 4 list of tetrahedra indices
  % Ouputs:
  %   v  #T list of tet volumes. Signed if dim = 3
  %

  a = V(T(:,1),:);
  b = V(T(:,2),:);
  c = V(T(:,3),:);
  d = V(T(:,4),:);
  % http://en.wikipedia.org/wiki/Tetrahedron#Volume
  % volume for each tetrahedron

  % Minus sign so that typical tetgen mesh has positive volume
  %v = -dot((a-d),cross2(b-d,c-d),2)./6./4;
  % Not sure where that ./4 came from...
  v = -dot((a-d),cross2(b-d,c-d),2)./6;
  function r = cross2(a,b)
    % Optimizes r = cross(a,b,2), that is it computes cross products per row
    % Faster than cross if I know that I'm calling it correctly
    r =[a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
        a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
        a(:,1).*b(:,2)-a(:,2).*b(:,1)];
  end
end
