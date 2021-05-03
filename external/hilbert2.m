function [x,y] = hilbert2(n)
% Hilbert 2D curve.
%
% function [x,y] = hilbert2(n) gives the vector coordinates of points
% in n-th order Hilbert curve of area 1.
%
% Example:
%   n=3;
%   [x,y] = hilbert2(n);
%   [V,E] = upsample([x(:) y(:)]*2^(n+1)+2^n,[1:numel(x)-1; 2:numel(x)]');
%   [sQ,sV] = surf2patch(zeros([1 1]*2^(n+1)));
%   I = accumarray(V,[(1:2:size(V,1)) 2:2:size(V,1)]',[1 1]*2^(n+1)-1);
%   sQ = sQ(logical(accumarray(V,true,[1 1]*2^(n+1)-1)),:);
%   sI = I(I>0);
%   tsurf(sQ,sV,'CData',sI);


% From hilbert3.m
%   Copyright (c) by Ivan Martynov
%   Inspired by function [x,y] = hilbert(n) made by Federico Forte
%   Date: September 17, 2009
% From hilber.m
%   Copyright (c) by Federico Forte
%   Date: 2000/10/06


if nargin ~= 1
    n = 2;
end

if n <= 0
    x = 0;
    y = 0;
else
    [xo,yo] = hilbert2(n-1);
    x=.5*[-.5+yo -.5+xo .5+xo  .5-yo];
    y=.5*[-.5+xo  .5+yo .5+yo -.5-xo];
end
