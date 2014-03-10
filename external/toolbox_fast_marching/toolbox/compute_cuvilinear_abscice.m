function s = compute_cuvilinear_abscice(c)

% compute_cuvilinear_abscice - compute the curvilinear abscice of a curve
%
%   s = compute_cuvilinear_abscice(c);
%
%   A curve is a ndims x npts matrix.
%
%   Copyright (c) 2004 Gabriel Peyré

if size(c,1)>size(c,2)
    c = c';
end

npts = size(c,2);
D = c(:,2:end)-c(:,1:(end-1));
s = zeros(npts,1);
s(2:end) = sqrt( sum(D.^2,1) );
s = cumsum(s);