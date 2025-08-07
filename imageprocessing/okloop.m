function [RGB,OK] = okloop(nc, arc, shift)
  % [RGB,OK] = okloop(nc, arc, shift)
  %
  % Example:
  %   r = 2; % number of repetitions
  %   n = 4; % number of unique colors
  %   okloop(r*n,r*2*pi);
  % 
  % Example:
  %   % n must be odd
  %   n = 15;
  %   X = okloop(n);
  %   X = X(mod((1:floor(n/2):n*floor(n/2)) -1,n)+1,:);
  %
  %   % similar hue range as jet
  %   okloop(256,-4/3*pi,-1/2*pi);
  if nargin < 1, nc = 256; end
  if nargin < 2, arc = 2*pi; end
  if nargin < 3, shift = 0; end

  th = linspace(0+shift,arc+shift,nc+1)';
  th = th(1:end-1);
  l = 0.75010101010101016;
  r = 0.12755316371916220;
  OK = [repmat(l,nc,1) r*[cos(th) sin(th)]];
  RGB = oklab2rgb(OK);
end
