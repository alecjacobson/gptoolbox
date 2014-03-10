function D = compute_levelset_shape(name, n, options)

% compute_levelset_shape - compute some basic level set shapes
%
%   D = compute_levelset_shape(name, n, options);
%
%   name can be: 'circle', 'rectangle', 'small-disks', 'circlerect1',
%   'circlerect2', 'square'.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;


[Y,X] = meshgrid(1:n,1:n);

switch lower(name)
    
    case 'circle'
        if isfield(options, 'center')
            center = options.center;
        else
            center = [n n]/2;
        end
        if isfield(options, 'radius')
            radius = options.radius;
        else
            radius = n/4;
        end
        D = sqrt( (X-center(1)).^2+(Y-center(2)).^2 ) - radius;
    case 'rectangle'
        if isfield(options, 'corner1')
            c1 = options.corner1;
        else
            c1 = [0.2 0.2]*n;
        end
        if isfield(options, 'corner2')
            c2 = options.corner2;
        else
            c2 = [0.8 0.8]*n;
        end
        D = max( max( c1(1)-X, X-c2(1) ), max( c1(2)-Y, Y-c2(2) ) );
    case 'circlerect1'
        if not(isfield(options, 'corner1'))
            options.corner1 = [0.2 0.2]*n;
        end
        if not(isfield(options, 'corner2'))
            options.corner2 = [0.8 0.8]*n;
        end
        options.center = options.corner2;
        options.radius = (options.corner2(1)-options.corner1(1))/2;
        D = max( -compute_levelset_shape('circle', n, options), ...
                 compute_levelset_shape('rectangle', n, options));
    case 'circlerect2'
        if not(isfield(options, 'corner1'))
            options.corner1 = [0.15 0.15]*n;
        end
        if not(isfield(options, 'corner2'))
            options.corner2 = [0.65 0.65]*n;
        end
        options.center = options.corner2;
        options.radius = (options.corner2(1)-options.corner1(1))/2;
        D = min( compute_levelset_shape('circle', n, options), ...
                 compute_levelset_shape('rectangle', n, options));
    case 'square'
        if isfield(options, 'width')
            width = options.width;
        else
            width = 0.8*n;
        end
        if isfield(options, 'center')
            center = options.center;
        else
            center = [1 1]*n/2;
        end
        options.corner1 = center(:)-repmat(width,2,1)/2;
        options.corner2 = center(:)+repmat(width,2,1)/2;
        D = compute_levelset_shape('rectangle', n, options);
        
    case 'small-disks'
        D = Inf*ones(n);
        if isfield(options, 'nbdisks')
            nbdisks = options.nbdisks;
        else
            nbdisks = 5;
        end
        x = linspace(0,n,nbdisks+1); x(end) = []; x = x+n/(2*nbdisks);
        options.radius = 0.75*n/(2*nbdisks);
        for i=1:nbdisks
            for j=1:nbdisks
                options.center = [x(i) x(j)];
                D = min(D, compute_levelset_shape('circle', n, options) );
            end
        end
    otherwise 
        error('Unknown shape.');
end
