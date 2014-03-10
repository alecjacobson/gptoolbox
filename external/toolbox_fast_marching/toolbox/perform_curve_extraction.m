function c_list = perform_curve_extraction(M, t, options)

% perform_curve_extraction - extract level set curves from an image.
%
%   c_list = perform_curve_extraction(M, t, options)
%
%   'M' is the image.
%   't' is the level.
%   'c_list' is a cell array of 2D curves.
%   'options' is an optional structure with fields: 
%       'max_nb' is the number of curves extracted (only the 'nb' longest).
%       'min_length' is the minimum length (in *pixels*) of extracted curves.
%
%   The image is always assumed to lie in [0,1]x[0,(ny-1)/(nx-1)]
%   where [nx,ny] = size(M).
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    t = 0;
end
if nargin<3
    options.null = 0;
end

[nx,ny] = size(M);
x = 0:1/(nx-1):1;
y = ( 0:1/(ny-1):1 )*( (ny-1)/(nx-1) );
c = contourc(x,y,M',[t,t]);
clear c_list;
c_list = {};
k = 0;
p = 1;
while p < size(c, 2)
    lc = c(2,p);   % length of the curve
    cc = c(:,(p+1):(p+lc));
    p = p+lc+1;
    k = k+1;
    c_list{ k } = cc;
end

n = length(c_list);

% filter by nb
if isfield(options, 'max_nb')
    max_nb = options.max_nb;
    l = zeros(n,1);
    for i=1:n
        l(i) = length(c_list{i});
    end
    [tmp,I] = sort(l);
    I = reverse( I );
    I = I( 1:(min(max_nb,n)) );
    c_list1 = c_list;
    c_list = {};
    for i = 1:length(I)
        c_list{i} = c_list1{I(i)};
    end        
end

% filter by size
if isfield(options, 'min_length')
    min_length = options.min_length;
    c_list1 = c_list;
    c_list = {};
    k = 0;
    for i=1:n
        l = length(c_list1{i});
        if l>=min_length
            k = k+1;
            c_list{k} = c_list1{i};
        end
    end
end
