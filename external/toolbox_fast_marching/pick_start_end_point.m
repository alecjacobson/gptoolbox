function [start_point,end_point,center_point] = pick_start_end_point(M, no_end_point, do_rounding)

% pick_start_end_point - pick start/end points for front propagation.
%
%   [start_point,end_point,center_point] = pick_start_end_point(M);
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<2
    no_end_point = 0;
end
if nargin<3
    do_rounding = 1;
end

clf;
hold on;
if ~isempty(M)
imagesc(M');
axis image;
axis off;
colormap gray(256);
end

m = size(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick center point
if nargout>2
    disp('Pick center point.');
    center_point = round( ginput(1) )';
    plot( [center_point(1),m(1)], [center_point(2) center_point(2)] );
    plot( center_point(1), center_point(2), 'ro' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick starting point
disp('Pick starting point.');
start_point = ginput(1)';
if do_rounding
    start_point = round(start_point);
end
if nargout>2
   start_point(2) = center_point(2);  % project 
end
plot( start_point(1), start_point(2), 'rx' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick ending point
if nargout>1
    if no_end_point
        end_point = start_point;
        end_point(2) = end_point(2)-3;
    else
        disp('Pick end point.');
        end_point = ginput(1)';
        if do_rounding
            end_point = round(end_point);
        end
        plot( end_point(1), end_point(2), 'rs' );
    end
end

hold off;