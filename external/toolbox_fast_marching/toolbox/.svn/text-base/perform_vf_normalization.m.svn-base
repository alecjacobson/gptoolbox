function v2 = perform_vf_normalization(v1)

% perform_vf_normalization - renormalize a vf.
%
%   v2 = perform_vf_normalization(v1);
%
%   Copyright (c) 2004 Gabriel Peyré


n = sum(v1.^2,3);
I = find(n<eps);
n(I) = 1;
v2 = prod_vf_sf(v1,1./sqrt(n));