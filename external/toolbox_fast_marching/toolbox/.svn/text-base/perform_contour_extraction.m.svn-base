function c = perform_contour_extraction(M, t)

% perform_contour_extraction - exctract the geometry of 
%   a black contour.
%
%   c = perform_contour_extraction(M, t);
%
%   the output is a sampled curve in [0,1]x[0,1]
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    t = mean(M(:));
end

if M(1)<t
    M = 2*t-M;
end

n = size(M,1);
x = 0:1/(n-1):1;

c = contourc(x,x,M,[t,t]);
lc = c(2,1);   % length of the curve
c = c(:,2:end);

% swap X/Y because curve is extracted in matrix coordinates
c = [c(2,:); c(1,:)];

if size(c,2)>lc+1
    warning('More than one curve detected.');
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old code

if nargin<2
    hv = 1;    
end

if hv == 2
    M = M';
end

N = size(M,1);
xcont = zeros(N,1);
ycont = zeros(N,1);

for i=1:N
    v = M(i,:);
    if v(1)==0
        I = find( v~=0 );
    else
        I = find( v==0 );
    end
    if length(I)==0
        xcont(i) = Inf;     % curve is undefined here
        ycont(i) = Inf;
    else
        xcont(i) = (i-1)/(N-1);
        ycont(i) = (I(1)-1)/(N-1);
    end
end


I = find(xcont~=Inf);
xcont = xcont(I);
ycont = ycont(I);

if nargin>2
    % use wavelet coef to realign the geometry
    if hv == 2
        MW = MW';
    end
    
    NW = size(MW,1);
    
    k = 0;
    t = 0;
    for i=1:NW
        v = MW(i,:);
        I = find(v==1);
        if length(I)>1
            xi = (i-1)/(NW-1);
            y1 = (I(1)-1)/(NW-1);
            y2 = (I(end)-1)/(NW-1);
            if y2-y1<0.2
                y = (y1+y2)/2;
                k = k+1;
                yi = interp1( xcont,ycont, xi, 'linear', 'extrap' ); 
                t = t + yi-y;
            end
        end
    end
    if k>0
        t = t/k;
        ycont = ycont - t;
    end
end


c = [xcont';ycont'];
