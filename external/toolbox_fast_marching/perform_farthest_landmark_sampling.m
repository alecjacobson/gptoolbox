function landmark = perform_farthest_landmark_sampling( W, landmark, DL_landmark, base_points, DL_base, nbr_points )

% perform_farthest_landmark_sampling - samples landmark using least error seeding strategy
%
% points = perform_farthest_landmark_sampling( W, points, nbr_points, DL, base_points );
%
%   points can be []
%   DL is the distance map to base_points, wich is a 2 x p matrix
%   
%   Copyright (c) 2005 Gabriel Peyré

if size(base_points,1)~=2
    base_points = base_points';
end
if size(landmark,1)~=2
    landmark = landmark';
end
if size(base_points,1)~=2
    error('base_points and landmark should be of size 2xp');
end

n = size(W,1);
nbr_bp = size(base_points,2);
nbr_lm = size(landmark,2);
% number of seeded point for each trial
nbr_trial = 50;


for i=1:nbr_points
    % new potential location at random
    trial = floor( rand(2,nbr_trial)*n ) + 1;
    err = [];
    for k=1:nbr_trial
        % compute distance to trial
        D = perform_fast_marching(W, trial(:,k));
        % new list of distance
        DL = cat(3, DL_landmark, D);
        % compute cummulative error with this new point distance 
        e = 0;
        for s=1:nbr_bp
            D1 = compute_distance_landmark(base_points(:,s),DL,landmark,'align');
            m = D1-DL_base(:,:,s);
            e = e + sum( m(:).^2 );
        end
        err = [err, e];
    end
    % retrieve best match
    [tmp,I] = min(err);
    % add it to list
    landmark = [ landmark, trial(:,I) ];
end