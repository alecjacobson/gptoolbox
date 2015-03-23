function interior_point = point_inside_polygon(V)
  % POINT_INSIDE_POLYGON  return some point inside of the given polygon
  %
  % interior_point = point_inside_polygon(V)
  %
  % Input
  %  V n by 2 vertex array
  % Output
  %  interior_point  a point inside V
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)

  assert(size(V,1) >= 3);
  % start with a point above the polygon in the middle
  p = 1.0 + [mean(V(:,1)) max(V(:,2))];

  % first hits are just mid points on segments
  first_hits = ...
    [V(:,1) + V([end, 1:(end-1)],1) , ...
    V(:,2) + V([end, 1:(end-1)],2)]/2;

  distance_to_p = ... 
    sum((repmat(p,size(first_hits,1),1) - first_hits).^2,2);

  first_hit = first_hits(distance_to_p == min(distance_to_p),:);
  first_hit = first_hit(1,:);

  x1 = p(1);
  y1 = p(2);
  x2 = first_hit(1);
  y2 = first_hit(2);
  x3 = V(1:end,1);
  y3 = V(1:end,2);
  x4 = V([end, 1:(end-1)],1);
  y4 = V([end, 1:(end-1)],2);
  % Find intersections along this ray with each line segment's line
  % from http://en.wikipedia.org/wiki/Line-line_intersection
  second_hits = ...
      [((x1.*y2-y1.*x2).*(x3-x4) - (x1-x2).*(x3.*y4-y3.*x4))./...
      ((x1-x2).*(y3-y4) - (y1 - y2).*(x3-x4)), ...
      ((x1.*y2-y1.*x2).*(y3-y4) - (y1-y2).*(x3.*y4-y3.*x4))./...
      ((x1-x2).*(y3-y4) - (y1 - y2).*(x3-x4))];

  % Check if intersection is really on line segment, handling vertical and
  % horizontal lines
  valid_second_hits = ...
    second_hits( ...
      ((second_hits(:,1) >= x3 & second_hits(:,1) <= x4 )| ...
      (second_hits(:,1) <= x3 & second_hits(:,1) >= x4 ) | ...
      (x3 == x4) ...
      ) & ...
      ((second_hits(:,2) >= y3 & second_hits(:,2) <= y4 )| ...
      (second_hits(:,2) <= y3 & second_hits(:,2) >= y4 ) | ...
      (y3 == y4) ...
      ) & ~((x3 == x4)&(y3 == y4))...
    ,:);

  distance_to_first_hit = ... 
    sum((repmat(first_hit,size(valid_second_hits,1),1) - ...
      valid_second_hits).^2,2);

  EPSILON = 10.0^-15;
  valid_second_hits = ...
    valid_second_hits(distance_to_first_hit>EPSILON,:);
  distance_to_first_hit = ...
    distance_to_first_hit(distance_to_first_hit>EPSILON,:);

  second_hit = valid_second_hits(...
    distance_to_first_hit == min(distance_to_first_hit),:);
  % shouldn't have to do this:
  second_hit = second_hit(1,:);

  interior_point = (first_hit + second_hit)/2;

end
