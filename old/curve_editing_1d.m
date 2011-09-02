function deformed = curve_editing_1d( fixed, undeformed, h)
  % Creates and solves system of equations to determine new positions of
  % points on a curve given target fixed positions
  %
  % INPUT
  % fixed       dict of fixed positions (index, value)
  % undeformed  original input positions
  % h           distance between elements
  %
  % OUTPUT
  % deformed    new positions given target fixed postions
  A, rhs = curve_editing_1d_system_fixed( fixed, undeformed, h);
  deformed = A \ rhs; 
end