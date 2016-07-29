function [P,E] = reuleaux_triangle(sam)
  % REULEAUX_TRIANGLE Construct a reuleaux triangle (curve of constant width)
  %
  % Inputs:
  %   sam  number of samples per side
  %%   th_max  portion of circle to use as side in radians (0,pi*2/3)
  % Outputs:
  %   P  #sam*3 by 2 list of vertex positions
  %   E  #E by 2 list of edge indices into P
  %
  if nargin<1
    sam = 20;
  end
  th_max = pi/3;
  th = linspace(0,th_max,sam+1)';
  E = [1:numel(th)-1;2:numel(th)]';
  P = [cos(th) sin(th)];
  %c = [0.5 sqrt(3)/6];
  c = 0.5*[P(1,:)+P(end,:)]+sqrt(3)/6*[P(end,:)-P(1,:)]*[0 1;-1 0];
  rot = @(ph) [cos(ph) -sin(ph);sin(ph) cos(ph)];
  R = rot(-2*pi/3);
  P = bsxfun(@minus,P(1:end-1,:),c);
  P = [P;P*R;P*R*R];
  P = P/sqrt(sum(P(1,:).^2));
  E = [1:size(P,1);2:size(P,1) 1]';
end
