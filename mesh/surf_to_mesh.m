function [V,F] = surf_to_mesh(X,Y,Z)
  % SURF_TO_MESH Convert a grid-based `surf` representation of surface to a
  % triangle mesh representation.
  %
  % [V,F] = surf_to_mesh(X,Y,Z)
  %
  % Inputs:
  %   X  m by n list of x-coordinates
  %   Y  m by n list of y-coordinates
  %   Z  m by n list of z-coordinates
  % Outputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indices
  %

  clean = true;
  s2p = surf2patch(X,Y,Z);
  F = s2p.faces;
  V = s2p.vertices;
  F = [F(:,[1 2 3]);F(:,[1 3 4])];

  if clean
    [V,VI,VJ] = remove_duplicate_vertices(V,1e-7);
    F = VJ(F);
    F = F(doublearea(V,F)>eps,:);
  end

end
