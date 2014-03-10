function [points,D] = perform_farthest_point_sampling( W, points, npoints, options )

% perform_farthest_point_sampling - samples points using farthest seeding strategy
%
% points = perform_farthest_point_sampling( W, points, npoints );
%
%   points can be [] or can be a (2,npts) matrix of already computed 
%       sampling locations.
%   
%   Copyright (c) 2005 Gabriel Peyre

options.null = 0;
if nargin<3
    npoints = 1;
end
if nargin<2
    points = [];
end
n = size(W,1);

aniso = 0;
d = nb_dims(W);
if d==4
    aniso = 1;
    d = 2; % tensor field
elseif d==5
    aniso = 1;
    d = 3; % tensor field
end
s = size(W);
s = s(1:d);

% domain constraints (for shape meshing)
L1 = getoptions(options, 'constraint_map', zeros(s) + Inf );
verb = getoptions(options, 'verb', 1);
mask = not(L1==-Inf);

if isempty(points)
    % initialize farthest points at random
    points = round(rand(d,1)*(n-1))+1;
    % replace by farthest point
    [points,L] = perform_farthest_point_sampling( W, points, 1 );
    Q = ones(size(W));
    points = points(:,end);
    npoints = npoints-1;
else
    % initial distance map
    [L,Q] = my_eval_distance(W, points, options);
%    L = min(zeros(s) + Inf, L1);
%    Q = zeros(s);
end

for i=1:npoints
    if npoints>5 && verb==1
        progressbar(i,npoints);
    end
    options.nb_iter_max = Inf;
    options.Tmax = Inf; % sum(size(W));
    %     [D,S] = perform_fast_marching(W, points, options);
    options.constraint_map = L;
    pts = points;
    if not(aniso)
       pts = pts(:,end);
    end
    D = my_eval_distance(W, pts, options);
    Dold = D;
    D = min(D,L); % known distance map to lanmarks
    L = min(D,L1); % cropp with other constraints
    if not(isempty(Q))
        % update Voronoi
        Q(Dold==D) = size(points,2);
    end
    % remove away data
    D(D==Inf) = 0;
    if isempty(Q)
        % compute farthest points
        [tmp,I] = max(D(:));
        [a,b,c] = ind2sub(size(W),I(1));
    else
        % compute farthest steiner point
        [pts,faces] = compute_saddle_points(Q,D,mask);
        a = pts(1,1); b = pts(2,1); c = [];
        if d==3
            c = pts(3,1);
        end
    end
    if d==2 % 2D
        points = [points,[a;b]];
    else    % 3D
        points = [points,[a;b;c]];
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D,Q] = my_eval_distance(W, x, options)

% D is distance
% Q is voronoi segmentation

options.null = 0;
n = size(W,1);
d = nb_dims(W);

if std(W(:))<eps
    % euclidean distance
    if size(x,2)>1
        D = zeros(n)+Inf;
        Q = zeros(n);
        for i=1:size(x,2)
            Dold = D; Qold = Q;
            D = min(Dold, my_eval_distance(W,x(:,i)));
            % update voronoi segmentation
            Q(:) = i;
            Q(D==Dold) = Qold(D==Dold);
        end
        return;
    end
    if d==2
        [Y,X] = meshgrid(1:n,1:n);
        D = 1/W(1) * sqrt( (X-x(1)).^2 + (Y-x(2)).^2 );
    else
        [X,Y,Z] = ndgrid(1:n,1:n,1:n);
        D = 1/W(1) * sqrt( (X-x(1)).^2 + (Y-x(2)).^2 + (Z-x(3)).^2 );
    end
    Q = D*0+1;
else
	[D,S,Q] = perform_fast_marching(W, x, options);
end