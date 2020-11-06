function [T] = cubic_cubic_intersect(CA,CB,tol)
  % CUBIC_CUBIC_INTERSECT Intersect two cubic Bezier curves.
  % 
  % [T] = cubic_cubic_intersect(CA,CB,tol)
  %
  % Inputs:
  %   CA  #4 by 2 list of control point locations
  %   CB  #4 by 2 list of control point locations
  % Outputs:
  %   T  #T by 2 list of parameters between [0,1] of intersects, so that
  %     intersect i occurs at T(i,1) on curve A and at T(i,2) on curve B.
  %
  % See also: cubic_eval, cubic_is_flat, cubic_split
  %

  if nargin<3
    tol = 1e-7;
  end
  
  %plot_cubic(CA);
  %hold on;
  %arrayfun(@(p) set(p,'Color','r'), plot_cubic(CB));
  %hold off;
  %axis equal;
  %axis([-1 1 -1 1]);
  %drawnow

  %if ~bounding_box_overlap(CA,CB)
  % Avoid function overhead
  if (any(max(CA)<min(CB)) || any(max(CB)<min(CA)))
    T = zeros(0,2);
    return;
  end
  A_flat = cubic_is_flat(CA,tol);
  B_flat = cubic_is_flat(CB,tol);
  % lazy just wait until they're both flat
  if A_flat && B_flat
    s = lineSegmentIntersect([CA(1,:) CA(4,:)],[CB(1,:) CB(4,:)]);
    if s.intAdjacencyMatrix
      T = [s.intNormalizedDistance1To2 s.intNormalizedDistance2To1];
    else
      T = zeros(0,2);
    end
    return;
  end

  [CA1,CA2] = cubic_split(CA,0.5);
  [CB1,CB2] = cubic_split(CB,0.5);
  T11 = cubic_cubic_intersect(CA1,CB1,tol)*0.5+0.5*[0 0];
  T12 = cubic_cubic_intersect(CA1,CB2,tol)*0.5+0.5*[0 1];
  T21 = cubic_cubic_intersect(CA2,CB1,tol)*0.5+0.5*[1 0];
  T22 = cubic_cubic_intersect(CA2,CB2,tol)*0.5+0.5*[1 1];
  T = [T11;T12;T21;T22];
end
