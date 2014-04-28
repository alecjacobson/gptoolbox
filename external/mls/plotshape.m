function plotshape(X, closed, symbol)
% PLOTSHAPE  This function allow to plot a single shape
%
% This function allow to plot a single shape represented in omogeneous
% coordinates or euclidean ones by a series of points
%
%  Params:
%
%  X		= [x; y] or [x; y; w] coord matrix for shape points
%  closed	= boolean that indicates if the shape must be closed
%  symbol	= symbol to be used in the plot
%
% Examples:
%
%  Plot a square
% >> plotshape([-1 -1;-1 1;1 1;1 -1]');

% Check args:
% Check number of arguments
if nargin<1
	error('Argument required for the points plotting');
end
if nargin<2
	closed = true;
end
if nargin<3
	symbol = '-';
end
% Check first argument
sizeX=size(X);
if sizeX(2)==0 return; end
if sizeX(1)<2 | sizeX(1)>3
	error('X must be a n*2 or n*3 matrix');
end

% Check the type of plot
pX=[];pY=[];
if sizeX(1)==2
	% Preparing the series of points
	for i=1:sizeX(2)
		pX = [pX X(1,i)];
		pY = [pY X(2,i)];
	end
else
	% Plot in omogeneous coords
	for i=1:sizeX(2)
		pX = [pX X(1,i)/X(3,i)];
		pY = [pY X(2,i)/X(3,i)];
	end
end

% Closing
if closed
	pX = [pX pX(1)];
	pY = [pY pY(1)];
end

% The real plot
plot(pX,pY,symbol);
