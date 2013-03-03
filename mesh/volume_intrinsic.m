function [vol] = volume_intrinsic(l)
  % VOLUME Compute volumes of tets defined intrinsically by edge lengths l
  %
  % v = volume(l)
  % 
  % Inputs:
  %   l  #T by 6 list of tetrahedra side lengths of edges opposite *face* pairs
  %     [23 31 12 41 42 43]
  % Ouputs:
  %   vol  #T list of tet volumes (always positive)
  %

  % http://en.wikipedia.org/wiki/Heron%27s_formula#Heron-type_formula_for_the_volume_of_a_tetrahedron

  % U, V, W, u, v, w are lengths of edges of the tetrahedron (first three form
  % a triangle; u opposite to U and so on)
  u = l(:,1); v = l(:,2); w = l(:,3); 
  U = l(:,4); V = l(:,5); W = l(:,6); 
  X = (w - U + v).*(U + v + w);
  x = (U - v + w).*(v - w + U);
  Y = (u - V + w).*(V + w + u);
  y = (V - w + u).*(w - u + V);
  Z = (v - W + u).*(W + u + v);
  z = (W - u + v).*(u - v + W);
  a = sqrt(x.*Y.*Z); 
  b = sqrt(y.*Z.*X); 
  c = sqrt(z.*X.*Y); 
  d = sqrt(x.*y.*z); 
  vol = sqrt( ...
    (-a + b + c + d).* ...
    ( a - b + c + d).* ...
    ( a + b - c + d).* ...
    ( a + b + c - d))./ ...
    (192.*u.*v.*w);

end

