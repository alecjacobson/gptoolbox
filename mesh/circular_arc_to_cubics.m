function [Pe,Ce] = circular_arc_to_cubics(Pnext,r,large_arc,sweep)
  % [Pe,Ce] = circular_arc_to_cubics(Pnext,r,large_arc,sweep)
  %
  % Inputs:
  %   Pnext  2d position of arc end
  %   r  radius in x direction
  %   large_arc  0 or 1
  %   sweep  0 or 1
  % Outputs:
  %   Pe  #Pe by 2 list of positions
  %   Ce  #Ce by 4 list of cubic bezier indices into rows of Pe
  %

  function [c,r] = circle_center(A, B, r, sweep)
    % Compute the midpoint M of A and B
    M = (A + B) / 2;
    % Compute the distance between A and M
    d = norm(A - M);
    % Calculate the distance from M to the center C
    if (r^2 - d^2)<0
      %warning('A,a radius too small');
      %r = d+eps(abs(d));
      r = d+1e-7;
    end
    h = sqrt(r^2 - d^2);
    % Determine the direction from A to B
    dir = (B - A) / norm(B - A);
    % Determine the direction perpendicular to AB
    perp_dir = [-dir(2), dir(1)];
    % Compute the two possible centers
    c = M + ((sweep==1)*1+(sweep~=1)*-1)  * h * perp_dir;
  end


  Pabs = [0 0];

  if isequal(Pabs,Pnext)
    Pe = zeros(0,2);
    Ce = zeros(0,4);
  end

  [c,r] = circle_center(Pabs, Pnext, r, xor(sweep,large_arc));
  [Pe,Ce] = ellipse_to_spline(c(1),c(2),r,r);

  is_flat_tol = 0.01;
  Pcut = [Pabs;Pnext];
  J = zeros(size(Pe,1),1);
  % this tol is in parametric distance
  snap_tol = 1e-7;
  for p = 1:size(Pcut,1)
    [~,Ie,Se,K] = point_spline_squared_distance(Pcut(p,:),Pe,Ce,is_flat_tol);
    if abs(Se - 0) < snap_tol
      %fprintf('snap %d to start: %g\n',p,Se);
      J(Ce(Ie,1)) = p;
    elseif abs(Se - 1)< snap_tol
      %fprintf('snap %d to end: %g\n',p,Se);
      J(Ce(Ie,4)) = p;
    else
      %fprintf('split %d at %g\n',p,Se);
      [Q1,Q2] = cubic_split(Pe(Ce(Ie,:),:),Se);
      J = [J(1:Ce(Ie,1));0;0;p;0;0;J(Ce(Ie,4):end)];
      Pe = [Pe(1:Ce(Ie,1),:);Q1(2:end,:);Q2(2:end-1,:);;Pe(Ce(Ie,4):end,:)];
      Ce = mod(((1:3:size(Pe,1)-1)'  + (0:3)) - 1,size(Pe,1))+1;
    end
  end
  
  if sum(J~=0) == 1
    % snapped to same point
    %just add a line, 0,⅓,⅔,1
    Pe = [Pabs; Pabs + (1/3)*(Pnext-Pabs); Pabs + (2/3)*(Pnext-Pabs); Pnext];
    Ce = [1 2 3 4];
    return;
  end


  % explicitly close up
  if isequal(Pe(1,:),Pe(end,:))
    Pe = Pe(1:end-1,:);
    J(1) = max(J(1),J(end));
    J(end) = [];
  end
  I1 = find(J==1);
  I2 = find(J==2);
  Ce = mod(((1:3:size(Pe,1)-1)'  + (0:3)) - 1,size(Pe,1))+1;

  l = edge_lengths(Pe,Ce(:,[1 4]));

  % Mid points
  BC = ...
    (0.5).^3.*Pe(Ce(:,1),:) + ...
    3.*(0.5).^3.*Pe(Ce(:,2),:) + ...
    3.*(0.5).^3.*Pe(Ce(:,3),:) + ...
    (0.5).^3.*Pe(Ce(:,4),:);
  larger = sum((BC-Pabs).*((Pnext-Pabs)*[0 1;-1 0]),2)>0;

  left = sum(l(larger));
  right = sum(l(~larger));
  if left < right
    larger = ~larger;
  elseif left == right
    % Not completely convinced. But this fragile adjustment may handle circles
    % made from two 180° arcs.
    larger = xor(larger,large_arc);
  end

  Ce = Ce(larger == large_arc,:);

  if ~any(Ce(:,1)==I1)
    Ce = fliplr(Ce);
  end
  assert( sum(Ce(:,1)==I1) == 1 )
  assert( sum(Ce(:,4)==I2) == 1 )
  E = [Ce(:,1:2);Ce(:,2:3);Ce(:,3:4)];
  A = sparse(E(:,1),E(:,2),1,size(Pe,1),size(Pe,1));
  I = [I1];
  while I(end) ~= I2
    I = [I;find(A(I(end),:))];
  end
  Pe = Pe(I,:);
  Ce = (1:3:size(Pe,1)-1)'  + (0:3);
end

