function O = orthocenter(V,F)
  % ORTHOCENTER Compute the orthocenter of a triangle mesh
  %
  % O = orthocenter(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  % Outputs:
  %   O  #F by 3 list of orthocenters
  %
  % See also: barycenter, incenter, circumradius
  %

  m = size(F,1);
  O = zeros(m,3);
  E = cat(3, ...
    V(F(:,3),:)-V(F(:,2),:), ...
    V(F(:,1),:)-V(F(:,3),:), ...
    V(F(:,2),:)-V(F(:,1),:));
  N = normals(V,F);
  P = cross(E,repmat(N,[1 1 3]),2);
  P = P./sqrt(sum(P.^2,2));

  [~,cosA] = internalangles(V,F);

  G = cat(3, ...
    cosA(:,2).*cosA(:,3), ...
    cosA(:,1).*cosA(:,3), ...
    cosA(:,1).*cosA(:,2));
  O = sum(P.*G,3);


end
