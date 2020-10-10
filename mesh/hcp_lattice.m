function C = hcp_lattice(I,J,K)
  % HCP_LATTICE generate the center locations of a hexagonal close-packed
  % lattice of unit-radius spheres.
  %
  % Inputs:
  %   I  #I list of x-coordinate indices
  %   J  #I list of y-coordinate indices
  %   K  #I list of z-coordinate indices
  % Outputs:
  %   C  #I by 3 list of center locations
  %
  % Example:
  %   [I,J,K] = meshgrid(0:9,0:9,0:9);
  %   C = hcp_lattice(I(:),J(:),K(:));
  %   [SV,SF] = lloyd_sphere(300);
  %   [V,F] = repmesh(SV,SF,C);
  i = I(:);
  j = J(:);
  k = K(:);
  % https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres#Simple_hcp_lattice
  C = [2*i + mod(j+k,2), sqrt(3)*(j+1/3*mod(k,2)), 2*sqrt(6)/3*k];
end
