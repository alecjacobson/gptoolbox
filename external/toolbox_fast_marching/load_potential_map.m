function [M,W] = load_potential_map(name, n, options)

% load_potential_map - prepare a potential map for FM computations
%
%   [M,W] = load_potential_map(name, n, options);
%
%   'M' is a cool 2D image you can use for display
%   'W' is the potential map
%
%   Copyright (c) 2004 Gabriel Peyré


options.null = 0;
file_name = ['data/', name, '.png'];

switch name
    case 'gaussian'
        x = -1:2/(n-1):1;
        [Y,X] = meshgrid(x,x);
        sigma = 0.4;
        M = exp( -(X.^2+Y.^2)/sigma^2 );
        M = 1-rescale( M );
        W = M + 0.01;
    case 'constant'
        M = ones(n)*0.5;
        W = M;
    case 'peaks'
        M = rescale( peaks(n) );
		W = M;
    case 'mountain';
        M = rescale( double( imread(file_name) ) );
        W = rescale( double(M) );
        W = W + 0.5;   
    case 'road2'
        M = rescale( double( imread(file_name) ) );
        M = rescale(double(M));
        if 0
            clf; imagesc(M); colormap gray(256);
            [x,y] = ginput(1); x = round(x); y = round(y);
            v = M(x,y);
        else
            v = 0.3059;
            v = 0;
        end
        W = abs( M-v );
        W = rescale(W)+0.001;  
    case 'peaks'
        M = rescale( peaks(n) );
    case 'cavern'
        M = rescale( double( imread(file_name) ) );
        W = M + 0.00001;
    case 'binary'
        M = ones(n); M(1:end/2,:) = .1;
        W = M;
    otherwise
        M = rescale( double( imread(file_name) ) );
        W = compute_edge_energy(M, 5, 0.05);
end

% remember the previous dynamic range
dr = [min(W(:)), max(W(:))];
% resize to get a square image
m = n/max(size(M));
m = round( size(M)*m );
warning off;
M = perform_image_resize(M,m);
W = perform_image_resize(W,m);
warning on;
W = rescale( W, dr(1), dr(2) );