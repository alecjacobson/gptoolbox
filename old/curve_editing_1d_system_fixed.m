function A,rhs = curve_editing_1d_system_fixed( fixed, undeformed, h)
  % Constructs the system matrix and right-hand side, whose solution is the
  % new positions of points on a curve given fixed target positions
  %
  % INPUT
  % fixed       dict of fixed positions (index, value)
  % undeformed  original input positions
  % h           distance between elements
  %
  % OUTPUT
  % A           system matrix
  % rhs         right-hand side
  
  A, rhs = curve_editing_1d_system(undeformed, h);
  
  % set the boundary conditions
  % dirichlet conditions for two points at each end
  rhs(1) = undeformed(1);
  rhs(2) = undeformed(2);
  rhs(end-1) = undeformed(end-1);
  rhs(end) = undeformed(end);
  
  % indentity-out corresponding rows in the system matrix
  I = eye(size(A));
  A(1,:) = I(1,:);
  A(2,:) = I(2,:);
  A(end-1,:) = I(end-1,:);
  A(end,:) = I(end,:);
  
  % set conditions for fixed points and indentity-out those rows
  rhs(fixed(1,:)) = fixed(2,:);
  A(fixed(1,:),:) = I(fixed(1,:));
  
end