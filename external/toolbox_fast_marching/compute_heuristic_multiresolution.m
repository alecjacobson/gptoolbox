function H = compute_heuristic_multiresolution( W, end_points, options )


if isfield(options, 'reduc_factor')
    reduc_factor = options.reduc_factor;
else
    reduc_factor = 0.3;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute coarse resolution map
if reduc_factor~=1
    Ws = perform_image_resize( W, size(W,1)*reduc_factor, size(W,2)*reduc_factor );
else
    Ws = W;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the heuristic map
clear options;
options.nb_iter_max = Inf;
start_points_small = round( end_points*reduc_factor );
if strcmp(propagation_type, 'normal')
    [Hs,S] = perform_fast_marching(Ws, start_points_small, options);
else
    center_point_small = round( center_point*reduc_factor );
    [Hs,S] = perform_circular_fast_marching_2d(Ws, start_points_small, center_point_small, options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extrapolate the heuristic
if reduc_factor~=1
    H = image_resize(Hs,size(W,1),size(W,2));
else
    H = Hs;
end