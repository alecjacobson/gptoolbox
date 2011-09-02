% [dx, dy] = imGradient(im, edges) 
%
% Compute the gradient of the image using smooth derivative filters
% optimized for accurate direction estimation.  Coordinate system
% corresponds to standard pixel indexing: X axis points rightward.  Y
% axis points downward.  EDGES specify boundary handling (see corrDn
% for options).
%
% Unlike matlab's new gradient function, which is based on local
% differences, this function computes derivatives using 5x5 filters
% designed to accurately reflect the local orientation content.

% EPS, 1997.
% original filters from Int'l Conf Image Processing, 1994.
% updated filters 10/2003: see Farid & Simoncelli, IEEE Trans Image Processing, 13(4):496-508, April 2004.
% Incorporated into matlabPyrTools 10/2004.

function [dx, dy] = imGradient(im, edges) 

if (exist('edges') ~= 1)
    edges = 'dont-compute';
end

%% kernels from Farid & Simoncelli, IEEE Trans Image Processing, 13(4):496-508, April 2004.
gp = [0.037659 0.249153 0.426375 0.249153 0.037659]';
gd = [-0.109604 -0.276691 0.000000 0.276691 0.109604]';

dx = corrDn(corrDn(im, gp, edges), gd', edges);
dy = corrDn(corrDn(im, gd, edges), gp', edges);

return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TEST:

%%Make a ramp with random slope and direction 
dir = 2*pi*rand - pi;
slope = 10*rand;

sz = 32
im = mkRamp(sz, dir, slope);
[dx,dy] = imGradient(im);
showIm(dx + sqrt(-1)*dy);

ctr = (sz*sz/2)+sz/2;
slopeEst = sqrt(dx(ctr).^2 + dy(ctr).^2);
dirEst = atan2(dy(ctr), dx(ctr));

[slope, slopeEst]
[dir, dirEst]
