function plotpoints(X,symbol)
% PLOTPOINTS  This function allow to plot a series of points
%
% This function allow to plot a series of points represented in omogeneous
% coordinates or euclidean ones
%
%  Params:
%
%  X		= [x; y] or [x; y; w] coord matrix (for multiple points)
%  symbol	= symbol to be used in the plot
%
% Warning:
%
%  This function use the hold on setting and mantain the old value unmodified
%
% Examples:
%
%  Plot of a pair of points as * (points [1;1] and [2;2])
% >> plotpoints([1 1; 2 2]','*')
% 
%  Plot the pair of points of the prevoius example expressed in omogeneous coords
% >> plotpoints([1 1 1; 20 20 10]');

% Check args:
% Check number of arguments
if nargin<1
	error('Argument required for the points plotting');
end
if nargin<2
	symbol = '+';
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
	% Plot in euclidean coords
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

% Plotting
plot(pX,pY,symbol);
