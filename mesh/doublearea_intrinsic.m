function dblA = doublearea_intrinsic(l)
  % DOUBLEAREA_INTRINSIC Compute the double area of the triangles of a mesh
  %
  % dblA = doublearea_intrinsic(l)
  %
  % Inputs:
  %   l  #F by 3, array of edge lengths of edges opposite each face in F
  % Outputs:
  %   dblA   #F list of twice the area of each corresponding face
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), and Daniele Panozzo
  %

  l = sort(l,2,'descend');
  l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);
  % Kahan's assertion: "Miscalculating Area and Angles of a Needle-like
  % Triangle" Section 2.
  % http://www.cs.berkeley.edu/~wkahan/Triangle.pdf
  if any(l3-(l1-l2)<0)
    warning( 'Failed Kahan''s assertion');
  end
  %% semiperimeters
  %s = (l1 + l2 + l3)*0.5;
  %% Heron's formula for area
  %dblA = 2*sqrt( s.*(s-l1).*(s-l2).*(s-l3));
  % Kahan's heron's formula
  dblA = 2*0.25*sqrt((l1+(l2+l3)).*(l3-(l1-l2)).*(l3+(l1-l2)).*(l1+(l2-l3)));
end


