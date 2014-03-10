function [D,S] = perform_circular_fast_marching_2d(W, start_points, center_point, options, H)

% perform_circular_fast_marching_2d - launch the Fast Marching algorithm for circular propagation.
%
%   [D,S] = perform_circular_fast_marching_2d(W, start_points,center_point,options, H);
%
%   'W' is the weight matrix (the highest, the slowest the front will move).
%   'start_points' is a 2 x k array, start_points(:,i) is the ith starting point .
%   'center_point' is a 2 x k array, it is the center around the flow will turn.
%
%   Optional:
%   - You can provide special conditions for stop in options :
%       'options.end_points' : stop when these points are reached
%       'options.nb_iter_max' : stop when a given number of iterations is
%          reached.
%   - You can provide an heuristic (typically that try to guess the distance
%       that remains from a given node to a given target).
%       This is an array of same size as W.
%
%   Copyright (c) 2004 Gabriel Peyré


if nargin<4
    options.null = 0;
end

if isfield(options,'end_points')
    end_points = options.end_points;
else
    end_points = [];
end

if isfield(options,'end_points')
    end_points = options.end_points;
else
    end_points = [];
end

if isfield(options,'verbose')
    verbose = options.verbose;
else
    verbose = 1;
end

if isfield(options,'nb_iter_max')
    nb_iter_max = options.nb_iter_max;
else
    nb_iter_max = Inf;
end

nb_iter_max = min(nb_iter_max, 1.2*max(size(W))^2);


% use fast C-coded version if possible
if exist('perform_front_propagation_2d')~=0
    if nargin<5 
        [D,S] = perform_circular_front_propagation_2d(W,start_points-1,end_points-1,center_point-1,nb_iter_max);
    else
        [D,S] = perform_circular_front_propagation_2d(W,start_points-1,end_points-1,center_point-1,nb_iter_max, H);
    end
else
    error('You must compile the mex files.');
end

% replace C 'Inf' value (1e9) by Matlab Inf value.
I = find( D>10000 );
D(I) = Inf;
