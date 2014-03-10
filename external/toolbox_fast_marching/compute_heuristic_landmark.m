function [H,Z] = compute_heuristic_landmark(DL,start_points)

% compute_heuristic_landmark - compute an heuristic using landmark points
%
%   [H,Z] = compute_heuristic_landmark(DL,start_points);
%
%   DL(:,:,i) is the distance map to the ith landmark point.
%   'H' is an approximation of the distance to 'start_point' (the heuristic map).
%
%   Copyright (c) 2005 Gabriel Peyré

n = size(DL,1);
T = repmat(DL(start_points(1),start_points(2),:), [n,n,1] );
[H,Z] = max( abs(T-DL), [], 3 );