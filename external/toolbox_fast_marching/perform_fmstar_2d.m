function [D,S] = perform_fmstar_2d(W, start_points,end_points, options)

% perform_fmstar_2d - launch the Fast Marching* algorithm.
%
%   [D,S] = perform_fmstar_2d(W, start_points,end_points, options)
%
%   'W' is the weight matrix (the highest, the slowest the front will move).
%   'start_points' is a 2 x k array, start_points(:,i) is the ith starting point .
%   'end_points' is a 2 x 1 array, it is the goal.
%
%   'options' is a structure that can contains the following fields: 
%       * 'heuristic': the value for the 'distance to goal' heuristic 
%           (should be of the same size as W)
%       * 'heuristic_type': if you do not provide explicitly the heuristic,
%           it can be computed automatically, and the way to compute
%           it are set by 'heuristic_type' (either 'multiresolution'
%           or 'landmarks')
%   'reduc_factor' is the reduction factor for the coarse scale
%       computation (eg. 0.5 will perform the heuristic computation
%       on a grid of half size).
%
%   Copyright (c) 2004 Gabriel Peyré

options.null = 0;
if isfield(options, 'weight')
    weight = options.weight;
else
    weight = 0.6;
end

if isfield(options, 'propagation_type')
    propagation_type = options.propagation_type;
else
    propagation_type = 'normal';
end
if strcmp(propagation_type, 'circular')
    if isfield(options, 'center_point')
        center_point = options.center_point;
    else
        warning('For circula mode, you must provide options.center_point: switching to normal mode.');
        propagation_type = 'normal';
    end
end

if isfield(options, 'heuristic')
    H = options.heuristic;
else
    if isfield(options, 'heuristic_type')
        heuristic_type = options.heuristic_type;
    else
        heuristic_type = 'multiresolution';
    end
    switch lower(heuristic_type)
        case 'multiresolution'
            H = compute_heuristic_multiresolution( W, end_points, options );
        case 'landmark'
            if isfield( options, 'distance_to_landmarks' )
                DL = options.distance_to_landmarks;
            elseif isfield( options, 'landmarks' )
                landmarks = options.landmarks;
                nbr_landmarks = size(nbr_landmarks, 2);
                % compute distance to landmark points
                DL = zeros(n,n,nbr_landmarks);
                opt.nb_iter_max = Inf; opt.Tmax = sum(size(W));
                for i=1:nbr_landmarks
                    [DL(:,:,i),S] = perform_fast_marching(W, landmarks(:,i), opt);
                end
            else
                error('For landmark heuristic, you should provide either the landmark points or the distance to them.');
            end
            H = compute_heuristic_landmark( DL, end_points );
        otherwise
            error('Unkwnown heuristic_type.');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the full FM
clear options;
options.nb_iter_max = Inf;
options.end_points = end_points;
if strcmp(propagation_type, 'normal')
    [D,S] = perform_fast_marching(W, start_points, options, H*weight);
else
    [D,S] = perform_circular_fast_marching_2d(W, start_points, center_point, options, H*weight);
end