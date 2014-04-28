function P = getpoints(N,force,mark)
% GETPOINTS  Allow to get from the current figure a set of N points
%
%  This function allow to select a set of N points by picing them with the
% mouse on the current figure. If the parameter is omitted the number of
% points is selected from the user.
%
%  Params:
%
% N         = The number of points. (def=any=0)
% force     = Points must be precisely N (true)? (def=true)
% mark      = The marker to be used, '' for no marker. (def='ro')
%
% P         = The points (as a matrix 2xN)
%
%  Examples:
%
%  Plot the selected points:
% >> imshow(img); hold on; plotpoints(getpoints);

% The def params:
if nargin<1 N = 0; end
if nargin<2 force = true; end
if nargin<3 mark = 'ro'; end

% Obtaining the points:
npts = 0;
P = [];
if N==0
        [X,Y] = getpts;
        P = [X,Y]';
        npts = numel(X);
else
    while (force && npts<N) || npts==0
        % Get another set of points:
        if strcmp('',mark)
            [X,Y] = ginput(N-npts);
        else
            [X,Y] = GetPoints(N-npts,mark);
        end

        % Converting the points format and counting:
        P = [P,[X,Y]'];
        npts = size(P,2);
    end
end

% If required removing the exceding points:
if nargin>0 && size(P,2)>N
    P = P(:,1:N);
end

% ------------------------ LOCAL FUNCTIONS ------------------------

% Get a set of points marking the points:
function [X,Y] = GetPoints(N,mark)

% Prepare holding
washold = ishold;
if ~washold
	hold on;
end

% The acquisition cicle:
X = zeros(N,1);
Y = zeros(N,1);
for ind=1:N
    % Get a point:
    [X(ind),Y(ind)] = ginput(1);
    % Mark the point:
    plot(X(1:ind),Y(1:ind),mark);
end

% Correct the old hold state
if ~washold hold off; end
