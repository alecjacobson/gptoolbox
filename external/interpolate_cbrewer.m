function [interp_cmap]=interpolate_cbrewer(cbrew_init, interp_method, ncolors)
% 
% INTERPOLATE_CBREWER - interpolate a colorbrewer map to ncolors levels
%
% INPUT:
%   - cbrew_init: the initial colormap with format N*3
%   - interp_method: interpolation method, which can be the following:
%               'nearest' - nearest neighbor interpolation
%               'linear'  - bilinear interpolation
%               'spline'  - spline interpolation
%               'cubic'   - bicubic interpolation as long as the data is
%                           uniformly spaced, otherwise the same as 'spline'
%   - ncolors=desired number of colors 
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 14.10.2011


% just to make sure, in case someone puts in a decimal
ncolors=round(ncolors);

% How many data points of the colormap available
nmax=size(cbrew_init,1);

% create the associated X axis (using round to get rid of decimals)
a=(ncolors-1)./(nmax-1);
X=round([0 a:a:(ncolors-1)]);
X2=0:ncolors-1;

z=interp1(X,cbrew_init(:,1),X2,interp_method);
z2=interp1(X,cbrew_init(:,2),X2,interp_method);
z3=interp1(X,cbrew_init(:,3),X2, interp_method);
interp_cmap=round([z' z2' z3']);

end
% Copyright (c) 2011, Charles Robert
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% 
