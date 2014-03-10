function M = perform_blurring(M, sigma, options)

% perform_blurring - gaussian blurs an image
%
%   M = perform_blurring(M, sigma, options);
%
%   M is the original data
%   sigma is the width of blurs (in pixels)
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
n = size(M,1);

eta = 4;
p = round((sigma*eta)/2)*2+1;
p = min(p,round(n/2)*2-1);

h = compute_gaussian_filter(p*[1 1],sigma/(2*n),[n n]);
M = perform_convolution(M, h, options);