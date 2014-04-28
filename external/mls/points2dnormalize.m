function Po = points2dnormalize(Pi,mustRemove)
% POINTS2DNORMALIZE  Normalize 2d points
%
%  Normalize 2d points removing the last vector element.
%
%  Params:
%
%  Pi           = Input points
%  mustRemove   = Must the third column be removed? (def=true)
%
%  Po           = Output points

% Check args:
% Check number of arguments
if nargin<1
	error('Argument required for the points to be normalized');
end

% Check the removal switch
if nargin<2
    mustRemove = true;
end

% Check the input points matrix geometry
sizePi=size(Pi);
if sizePi(2)==0 Po=Pi;return; end
if sizePi(1)<2 | sizePi(1)>3
	error('Pi must be a n*2 or n*3 matrix');
end

% Check for normalization:
Po = Pi;
if not(sizePi(1)==2)
    % Normalization
    for i=1:sizePi(2)
        Po(:,i) = Po(:,i)/Po(3,i);
    end
end

% Returning
if mustRemove
    Po = Po(1:2,:);
elseif sizePi(1)==2
    Po = [Pi;ones(1,size(Pi,2))];
end
