function W = compute_edge_energy(M,s,epsi, center_point)

% compute_edge_energy - compute an energy for fast marching.
%
%   W = compute_edge_energy(M,s,epsi);
%
%   W is the speed function for the front propagation
%   (should be high in the area of strong gradient).
%   The formula is :
%
%   1/W(x) = 1/(1 + |grad(M*G_s)|) + epsi
%
%   where G_s is a gaussian smoothing of width s.
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<3
    epsi = 0.05;
end
if nargin<2
    s = 3;
end

M = rescale(M);
s = round(s);

if ndims(M)==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % smooth a bit
%    h = ones(s)/s^2;
    h = compute_gaussian_filter([15 15], s/(2*size(M,1)), size(M));
    G = conv2(M,h, 'same');
    grad = compute_grad(G);
    G = 1./(1+sqrt(sum(grad.^2,3))) + epsi;
    W = 1./G;        % speed = 1/potential
    % remove border artifacts
    e = 1/(1+epsi);
    W(1:s+2,:) = e;
    W(:,1:s+2) = e;
    W(end-s-1:end,:) = e;
    W(:,end-s-1:end) = e;

    if nargin>3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % scale the metric to avoid shrinkage near zeros
        m = size(M);
        [Y,X] = meshgrid(1:m(2), 1:m(1));
        d = sqrt( (X-center_point(1)).^2 + (Y-center_point(2)).^2 );
        d(center_point(1),center_point(2)) = 1e-9;
        W = W.*d;
    end
elseif ndims(M)==3
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [a,b,c] = size(M);
    % smooth a bit
    G = smooth3(M,'box',s);
    [gx,gy,gz] = gradient(G,1/a,1/b,1/c);
    G = 1./(1+sqrt(gx.^2 + gy.^2 + gz.^2)) + epsi;
    W = 1./G;        % speed = 1/potential

    if nargin>3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % scale the metric to avoid shrinkage near zeros
        m = size(M);
        [X,Y,Z] = meshgrid(1:m(2), 1:m(1));
        d = sqrt( (X-center_point(1)).^2 + (Y-center_point(2)).^2 + (Y-center_point(3)).^2 );
        d(center_point(1),center_point(2),center_point(3)) = 1e-9;
        W = W.*d;
    end
else
    error('Works only for 2D and 3D data.');
end