function plotpointsLabeled(X,symbol,labels,color)
% PLOTPOINTSLABELED  Plot a set of points labeling them
%
%  This function allow to plot a series of points represented in omogeneous
% coordinates or euclidean ones labeling with numeric or arbitrarious text
% the points.
%
%  Params:
%
%  X		= [x; y] or [x; y; w] coord matrix (for multiple points)
%  symbol	= Symbol to be used in the plot
%  labels   = Labels used for the points, must be a string cell of N
%             strings, starting number is accepted. (def={}=numbers)
%  color    =  The color used for the labels as [r,g,b] with range [0:1]
%             for each channel. (def=[0,0,0])
%
% Examples:
%
%  Plot a set of numbered points
% >> plotpointsLabeled([1 1; 2 2]','*')

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
% Check the labels:
if nargin<3 || ~isa(labels,'cell') || numel(labels)==0
    % Init the generation:
    if nargin>=3 && ~isa(labels,'cell') off=labels;
    else off=1; end
    labels = {};
    % Generating numerical labels:
    for ind=1:sizeX(2)
        labels{ind} = sprintf('%d',ind+off-1);
    end
end
% The color:
if nargin<4 color=[0,0,0]; end

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

% Plotting labels
for i=1:min(sizeX(2),size(labels,2))
    % Plotting a single label
    text(pX(i),pY(i),labels(i),'Color',color);
end
