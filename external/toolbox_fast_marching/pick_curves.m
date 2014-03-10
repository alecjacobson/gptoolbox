function sk = pick_curves(mask)

% pick_curves - ask for the user to build a set of curves
%
%   sk = pick_curves(mask);
%
%   mask is a background image (should be in [0,1] approx).
%
%   The user right-click on a set of point which create a curve.
%   Left click stop a curve.
%   Another left click stop the process.
%
%   Copyright (c) 2007 Gabriel Peyr?


n = size(mask,1);

sk = zeros(n);

b = 1;
while b(end)==1
    % draw a line
    clf;
    imagesc(mask+sk); axis image; axis off;
    colormap gray(256);
    [y1,x1,b] = ginput(1);
    while b==1
        clf;
        imagesc(mask+sk); axis image; axis off;
        [y2,x2,c] = ginput(1);
        if c~=1
            break;
        end
        t = linspace(0,1,1000);
        x = round( x1*t+x2*(1-t) );
        y = round( y1*t+y2*(1-t) );
        sk( x+(y-1)*n ) = 1;
        x1 = x2; y1 = y2;
    end
end