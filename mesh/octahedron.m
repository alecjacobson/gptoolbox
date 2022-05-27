function [V,F] = octahedron()
  V = [eye(3);-eye(3)];
  F = [
     1     2     3
     1     3     5
     1     5     6
     1     6     2
     2     4     3
     2     6     4
     3     4     5
     4     6     5];
end
