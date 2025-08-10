function [V,E] = isocontour(X,Y,Z,v)
  % [V,E] = isocontour(X,Y,Z,v)
  %
  % Wrapping on contourc to return isocontour as edges and vertices.
  %
  % Inputs:
  %   X  ny by nx matrix of x-coordinates
  %   Y  ny by nx matrix of y-coordinates
  %   Z  ny by nx matrix of function values
  %   v  scalar value for the isocontour
  % Outputs:
  %   V  #V by 2 list of vertex positions
  %   E  E# by 2 list of edges as indices into rows of V
  c = contourc(X(1,:),Y(:,1),Z,v);
  i = 1;
  V = [];
  E = [];
  I = [];
  L = [];

  while i<=size(c,2)
    if isequal(c(:,i+1),c(:,i+c(2,i)))
      if c(2,i) > 1
        E = [E size(V,2)+ [1:c(2,i)-1;2:c(2,i)-1 1]];
        V = [V c(:,i+1:i+c(2,i)-1)];
        I = [I repmat(numel(L),1,c(2,i)-1)];
      end
    else
      E = [E size(V,2)+ [1:c(2,i)-1;2:c(2,i)]];
      V = [V c(:,i+1:i+c(2,i))];
      I = [I repmat(numel(L),1,c(2,i))];
    end
    L = [L c(1,i)];
    i = i+1+c(2,i);
  end
  V = V';
  E = E';
end

