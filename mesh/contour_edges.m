function [V,E] = contour_edges(X,Y,Z,v)
  % CONTOUR_EDGES Extract edges of a contour from a 2D grid
  %
  % [V,E] = contour_edges(X,Y,Z,v)
  %
  % Input:
  %   X  h by w matrix of x coordinates
  %   Y  h by w matrix of y coordinates
  %   Z  h by w matrix of values
  %   v  value to extract
  % Output:
  %   V  n by 2 matrix of vertices
  %   E  m by 2 matrix of edges
  %

  %% This doesn't orient correctly
  %[cc] = contourc(X(1,:),Y(:,1),Z,[v v]);
  %pointer = 1;
  %V = zeros(2,0);
  %E = zeros(2,0);
  %while true
  %  if pointer > size(cc,2)
  %    break;
  %  end
  %  n = cc(2,pointer);
  %  x = cc(1,pointer+1:pointer+n);
  %  y = cc(2,pointer+1:pointer+n);
  %  if x
  %    E = [E size(V,2)+[1:n-1;2:n]];
  %  else
  %    n = n-1;
  %    x = x(1:n);
  %    y = y(1:n);
  %    E = [E size(V,2)+[1:n;2:n 1]];
  %  end
  %  V = [V [x;y]];

  %  pointer = pointer + n + 1;
  %end
  %V = V';
  %E = E';

  [E,V] = isocontour_off(Z,v);
  E = fliplr(E);

  V = V(:,[2 1]);
  V = V - 1;
  dims = size(X);
  XY = [X(:) Y(:)];
  V = V./(dims([2 1])-1).*(max(XY)-min(XY))+min(XY);


end
