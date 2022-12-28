function q = super_fibonacci(n,phi,psi)
  % SUPER_FIBONACCI Generate n quaternions according to "Super-Fibonacci
  % Spirals: Fast, Low-Discrepancy Sampling of SO(3)" [Alexa 2021]
  %
  % q = super_fibonacci(n,phi,psi)
  % 
  % Inputs:
  %   n  number of rotations to generate
  %   phi, psi  parameters controling the distribution {sqrt(2), 1.53375116875…}
  % Outputs:
  %   q  n by 4 list of unit quaternions
  %  
  % 
  if nargin<2
    phi = sqrt(2);
  end
  if nargin<3
    psi = 1.533751168755204288118041;
  end
  % https://marcalexa.github.io/superfibonacci/
  i = (0:n-1)';
  s = i + 0.5;
  r = sqrt(s./n);
  R = sqrt(1.0-s./n)
  alpha = (2*pi .* s)./phi;
  beta = (2*pi .* s)./psi;
  q = [r.*sin(alpha), r.*cos(alpha), R.*sin(beta), R.*cos(beta)];
end
