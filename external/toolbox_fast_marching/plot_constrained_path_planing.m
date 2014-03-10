function plot_constrained_path_planing(M,path,obj,nb_s,options)

% plot_constrained_path_planing - plot the result of the fast marching.
%
%   plot_constrained_path_planing(M,path,obj);
%
%   Copyright (c) 2004 Gabriel Peyré

options.null = 1;
if nargin<4
    nb_s = 20;
end
if isfield(options, 'box_width')
    box_width = options.box_width;
else
    box_width = 2;
end

n = size(M,1);


hold on;
imagesc(1-M');
colormap gray;
axis image;
axis square;
set(gca,'box','on');
if nargin>=3
%    plot(path(:,1), path(:,2), '-');
end

% resample the path evenly
path1 = perform_curve_resampling(path, nb_s, 'nbpts')';

% plot tiny red boxes
for i=1:size(path1,1)
    plot_box( path1(i,1:2), obj*n, path1(i,3), 'r', box_width );
end

set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca,'box','on');
hold off;

function plot_box(c,obj,theta, str, box_width)

v = [cos(theta),sin(theta)]*obj(1)/2;
w = [-sin(theta),cos(theta)]*obj(2)/2;
plot_line( c+v+w, c-v+w, str, box_width );
plot_line( c-v+w, c-v-w, str, box_width );
plot_line( c-v-w, c+v-w, str, box_width );
plot_line( c+v-w, c+v+w, str, box_width );

function plot_line(a,b, str, box_width)

h = line([a(1) b(1)], [a(2) b(2)], 'LineWidth', box_width);
set(h,'Color','r');