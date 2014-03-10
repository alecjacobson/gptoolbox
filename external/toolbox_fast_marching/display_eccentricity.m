function display_eccentricity(E, mode)

% display_eccentricity - display the transform
%
% display_eccentricity(E, mode);
%
%   mode==1 for colors
%   mode==2 for level sets
%
%   Copyright (c) 2006 Gabriel Peyré

if nargin<2
    mode = 1;
end


if mode==1
    I = find(E>0); A = E*0+256;
    A(I) = rescale(E(I),1,255);
    image(A);
    a = jet(256);
    a(end,:) = [1 1 1];
    colormap(a);
    axis image; axis off;
else
    s = unique(E(:)); I = find(s>0.1); s = s(I);
    x = linspace(s(1), s(end), 40);
    contour(E, x);
    axis image; axis off; axis ij;
end