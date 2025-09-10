function [V,E,F] = circle(n)
  th = linspace(0, 2*pi, n+1)';
  th = th(1:end-1); % Remove the last point to avoid duplication
  V = [cos(th), sin(th)]; % Circle in the XY plane
  E = [1:n; 2:n 1]'; % Edges connecting points in a circular manner
  F = [2:n-1; 3:n; ones(1,n-2)]'; % Faces defined by the edgesh
end
