function [V,F,UV] = moebius_strip(u,v)
  % MOEBIUS_STRIP  Generate a moebius strip
  %
  % Inputs:
  %   u  Number of vertices in u direction
  %   v  Number of vertices in v direction
  % Outputs:
  %   V  u*v by 3 list of vertex positions
  %   F  #F by 3 list of faces
  %   UV u*v by 3 list of original UV coordinates
  %   


  [UV,F] = create_regular_grid(u+1,v,0,0);
  % Reparameterize along x
  V = [(1-cos(UV(:,1)*pi))/2 UV(:,2)];
  % Twist once and smooth height
  V = [ ...
    1-acos(V(:,1)*2-1)/pi ...
    V(:,2).*V(:,1)+(1-V(:,2)).*(1-V(:,1)) ...
    (1-cos(V(:,1)*2*pi))/2.*(0.5-V(:,2))];
  % Wrap around
  V = [cos(V(:,1)*2*pi).*(2+V(:,2)) sin(V(:,1)*2*pi).*(2+V(:,2)) V(:,3)];
  I = 1:size(V,1);
  I(u*v + (1:v)) = v:-1:1;
  F = I(F);
  V = V(1:u*v,:);
  UV = UV(1:u*v,:);
end
