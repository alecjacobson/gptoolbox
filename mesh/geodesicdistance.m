function [ d ] = geodesicdistance( V,F, SP)
% [d] = geodesicdistance(V,F,SP)
%
% Compute the Geodesic Distance using the toolbox_fast_marching
% http://www.ceremade.dauphine.fr/~peyre/
%
% Input:
% V,F: mesh definition
% SP: an array of starting points
%
% Output:
% d: geodesic distance from the closest seed point computed on vertices,
% Vx1

options.end_points = [];
options.nb_iter_max = Inf;
d = perform_fast_marching_mesh(V, F, SP, options);
end

