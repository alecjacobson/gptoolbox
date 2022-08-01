function Q = fibonacci_sphere_sampling(P,phi)
  % FIBONACCI_SPHERE_SAMPLING Generate P points on the unit sphere using a
  % fibonacci spiral.
  %
  % Q = fibonacci_sphere_sampling(P,phi)
  % 
  % Inputs:
  %   P  number of points to generate
  %   phi  golden ratio {(1+sqrt(5))/2}
  % Outputs:
  %   Q  P by 3 list of points on the unit sphere
  %   
  if nargin<2
    phi = (1+sqrt(5))/2;
  end
  % https://arxiv.org/pdf/0912.4540.pdf
  % P = 2*N+1
  N = (P-1)/2;
  i = linspace(-N,N,P)';
  sin_lat = 2.*i./P;
  cos_lat = sqrt(1-sin_lat.^2);
  lon = 2.*pi.*i./phi^2;
  Q = [cos_lat.*cos(lon) cos_lat.*sin(lon) sin_lat];
end
