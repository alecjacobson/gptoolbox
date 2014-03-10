function [D,S] = perform_fmstar_3d(W, start_points,end_points, options)

% perform_fmstar_3d - launch the Fast Marching* algorithm in 3D.
%
%   [D,S] = perform_fmstar_2d(W, start_points,end_points, options)
%
%   'W' is the 3D weight matrix (the highest, the slowest the front will move).
%   'start_points' is a 3 x k array, start_points(:,i) is the ith starting point .
%   'end_points' is a 3 x 1 array, it is the goal.
%
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
if isfield(options, 'reduc_factor')
    reduc_factor = options.reduc_factor;
else
    reduc_factor = 0.3;
end

if reduc_factor~=1
    Ws = image_3d_resize(W,round(size(W)*reduc_factor));
else
    Ws = W;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the heuristic map
clear options;
options.nb_iter_max = 10000000;
start_points_small = round( end_points*reduc_factor );
[Hs,S] = perform_fast_marching(Ws, start_points_small, options);
% extrapolate the heuristic
if reduc_factor~=1
    H = image_3d_resize(Hs,size(W));
else
    H = Hs;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the full FM
clear options;
options.nb_iter_max = 10000000;
options.end_points = end_points;
[D,S] = perform_fast_marching(W, start_points, options, H*weight);