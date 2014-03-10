function M1 = generate_constrained_map(M,nb_theta,obj,theta_max)

% generate_constrained_map - generate a constraint map from a 2D map.
%
%   M1 = generate_constrained_map(M,nb_theta,obj);
%
%   M is the 2D map (binary image, 0 for background, 1 for objects).
%   M1 is a 3D array, M1(:,:,k) is the 2D map for orientation
%       theta=(t-1)*2*pi/nb_theta.
%  obj is the size of the size of the rectangle object that moves
%   (eg. obj=[0.2,0.05]).
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    nb_theta = 10;
end

if nargin<3
    obj = [0.2,0.05];
end

if nargin<4
    theta_max = 2*pi;
end

[n,p] = size(M);
M1 = zeros(n,p,nb_theta);

w = max(obj);
nw = norm(obj)*n+4;
nw = round((nw-1)/2)*2+1;

x = (-(nw-1)/2:(nw-1)/2)*1/(n-1);
[Y,X] = meshgrid(x,x);
for t = 1:nb_theta
    theta = theta_max*(t-1)/nb_theta;
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    Xr = X*cos(theta) + Y*sin(theta);
    Yr =-X*sin(theta) + Y*cos(theta);
    SE = (abs(Xr)<obj(1)/2) .* (abs(Yr)<obj(2)/2);
    M1(:,:,t) = imdilate(M,SE);
end