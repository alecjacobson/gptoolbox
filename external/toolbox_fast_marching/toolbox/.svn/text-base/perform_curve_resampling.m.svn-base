function c1 = perform_curve_resampling(c, s, mode, interp_type)

% perform_curve_resampling - resample a curve
%   
%   c1 = perform_curve_resampling(c, n, 'nbpts', intep_type);
% or
%   c1 = perform_curve_resampling(c, s, 'dist', intep_type);
%
% Where 'n' is the new number of point in the curve, 
% and 's' is the curvilinear absice between two points.
%
% For 'interp_type', see the 'method' parameter in 'help interp1'.
%
%   Copyright (c) 2005 Gabriel Peyré


if nargin<3
    mode = 'nbpts';
end
if nargin<4
    interp_type = 'linear';
end
    

if size(c,2)<size(c,1)
    c = c';
end


cabs = compute_cuvilinear_abscice(c);

if strcmp(mode, 'dist')
    n = round(cabs(end)/s)+1;
else
    n = s;
end

cabs_i = ( 0:cabs(end)/(n-1):cabs(end) )';



c1 = zeros(size(c,1), n);
for k=1:size(c,1)
    c1(k,:) = interp1( cabs, c(k,:), cabs_i );
end