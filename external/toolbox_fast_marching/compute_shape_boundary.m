function bound = compute_shape_boundary(M)

% compute_shape_boundary - extract boundary points
%
% bound = compute_shape_boundary(M);
%
%   bound is oriented counter clockwise.
%   This is the 8 connectivity boundary.
%
%   Copyright (c) 2007 Gabriel Peyre


n = size(M,1);

% enforce only 1 connected component
epsilon = 1e-10;
options.nb_iter_max = Inf;
options.Tmax = 2*n;
W = epsilon+M;
I = find(M(:)==1); x = linspace(-1,1,n);
[Y,X] = meshgrid(x,x); 
[tmp,c] = min( X(I).^2 + Y(I).^2 );
[a b] = ind2sub([n n],I(c));
options.constraint_map = [];
[d,tmp] = perform_fast_marching(W, [a;b], options);
I = find(d>1e5 | d<0); M(I) = 0;

% points on the boundary
if 0
    % outside the shape
    h = [0 1 0; 1 1 1; 0 1 0];
    Mh = perform_convolution(M,h);
    A = Mh>0 & M==0;
else
    % inside the shape
    h = ones(3);
    Mh = perform_convolution(M,h);
    A = Mh>0 & Mh<9 & M==1;
end

% 8 connectivity
wx = [1 -1 0  0 -1 1 -1  1];
wy = [0  0 1 -1 -1 1  1 -1];


% order the points
a = find(A(:)==1); a = a(1);
[x,y] = ind2sub(size(A),a);
    
boundx = x;
boundy = y;
A(x,y) = 0;


while sum(A(:))>0
    % find closest point
    for k=1:8
        if x+wx(k)>0 && x+wx(k)<=n && y+wy(k)>0 && y+wy(k)<=n && A(x+wx(k),y+wy(k))>0
            x = x+wx(k); 
            y = y+wy(k);
            break;
        end
    end
    if A(x,y)==0
        % try to find the existing point the closest to the last one
        I = find(A>0);
        [xlist,ylist] = ind2sub(size(A),I);
        d = (xlist-x).^2 + (ylist-y).^2;
        [tmp,I] = min(d(:));
        x = xlist(I(1)); y = ylist(I(1));
    end
    A(x,y) = A(x,y)-1;
    boundx(end+1) = x;
    boundy(end+1) = y;
end
bound = cat(2,boundx(:), boundy(:));


% reorient the curve so that it is clockwise oriented
p = size(bound,1);
c = mean(bound,1);
v = bound - repmat(c, [p 1]);
% cross product
z = v(1:end-1,1).*v(2:end,2) - v(1:end-1,2).*v(2:end,1);
if sum(z<0)>sum(z<0)
    % should flip the curve
    bound = bound(end:-1:1,:);
end