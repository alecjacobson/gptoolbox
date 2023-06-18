function [Pe,Ce] = arc_to_cubics(Pabs,Pnext,rx,ry,phi,large_arc,sweep)
  % [Pe,Ce] = arc_to_cubics(Pabs,Pnext,rx,ry,phi,large_arc,sweep)
  %
  % Inputs:
  %   Pabs  2d position of arc start
  %   Pnext  2d position of arc end
  %   rx  radius in x direction
  %   ry  radius in y direction
  %   phi  rotation of x axis in degrees
  %   large_arc  0 or 1
  %   sweep  0 or 1
  % Outputs:
  %   Pe  #Pe by 2 list of positions
  %   Ce  #Ce by 4 list of cubic bezier indices into rows of Pe
  %

  if rx == 0 || ry == 0
    warning('arc_to_cubics:zero-radius rx or ry is zero');
    Pe = zeros(0,2);
    Ce = zeros(0,4);
    return;
  end

  phi_rad = phi*pi/180;
  R = [cos(phi_rad) -sin(phi_rad) 0;sin(phi_rad) cos(phi_rad) 0;0 0 1];
  S = diag([1 rx/ry 1]);
  T = [eye(2,3);-Pabs 1] * R * S;

  [Pe,Ce] = circular_arc_to_cubics([Pnext 1]*T(:,1:2),rx,large_arc,sweep);
  Pe = [Pe ones(size(Pe,1),1)] * inv(T);
  Pe = Pe(:,1:2);

end
