function [flag, t, lambda] = ray_polygon_intersect(o,d,V,E)
% RAY_POLYGON_INTERSECT 2D Ray/polygon intersection
% 
% [flag, t, lambda] = ray_polygon_intersect(o,d,V,E)
%
% Input:
%    o  2D vector ray origin.
%    d  2D vector ray direction.
%    V  #V by 2 list of vertex positions
%    E  #E by 2 list of edge indices
% Output:
%    flag  #E list of bools: (false) Reject, (true) Intersect.
%    t  #E list of distances from the ray origin.
%    lambda  #E list of parameter of hit between E(:,1) and E(:,2)
%

  epsilon = eps;%0.00001;

  assert(size(V,2) == 2);
  assert(size(E,2) == 2);
  % number of edges 
  m = size(E,1);
  assert(numel(d) == 2);
  % make direction a row vector
  d = reshape(d,1,2);
  %d = normalizerow(d);
  d = d ./ sqrt(sum(d.^2,2));
  assert(numel(o) == 2);
  % make origin a row vector
  o = reshape(o,1,2);

  p1 = V(E(:,1),:);
  p2 = V(E(:,2),:);
  % edge vectors from 1 to 2
  e12 = p2-p1;
  % perpendiculars of edge vectors
  pe12 = perp(e12);


  % p1 minus o
  p1mo = plusrow(p1,-o);
  % pe12 dot d
  pe12dd = dotrow(pe12,d);
  % perp d
  pd = perp(d);

  % project to edges
  % http://objectmix.com/graphics/132701-ray-line-segment-intersection-2d.html
  t = dot2(pe12,p1mo)./pe12dd;
  lambda = dotrow(plusrow(-p1,o+d),pd)./dotrow(e12,pd);
  %lambda = dotrow(perp(d),p1mo)./pe12dd;

  flag = true(m,1);
  flag(lambda > 1) = 0;
  flag(lambda < 0) = 0;
  flag(t < 0) = 0;
  
  function u = perp(v)
    % [x,y] = [y,-x]
    u = [v(:,2),-v(:,1)];
  end

  function r = dot2(a,b)
    % Optimizes r = dot(a,b,2), that is it computes dot products per row
    % Faster than dot if I know that I'm calling it correctly
    r = sum(a.*b,2);
  end

  function r = dotrow(a,b)
    % Computes dot product of rows in a against b.
    % optimizes r = dot(a,repmat(b,size(a,1),1),2)
    r = a(:,1).*b(1,1) + a(:,2).*b(1,2);
  end

  function r = plusrow(a,b)
    % Computes sum rows in a with b.
    % optimizes r = a + repmat(b,size(a,1),1)
    r = [a(:,1)+b(1,1) a(:,2)+b(1,2)];
  end

end
