function plot_volumetric_data(W, nb_colors)

% plot_volumetric_data - plot a cube of data.
%
%   plot_volumetric_data(W, nb_colors);
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<2
    nb_colors = 256;
end

[n,p,q] = size(W);
h = slice(W,p/2,n/2,q/2);
set(h,'FaceColor','interp','EdgeColor','none');

box on;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'ZTick', []);

colormap( jet(nb_colors) );

lighting phong;
camlight infinite; 
camproj('perspective');


view(3);
axis tight;
axis equal;
% cameramenu;