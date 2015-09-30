function [colormap]=cbrewer(cname, ncol, interp_method)
% CBREWER - This function produces a colorbrewer table (rgb data) for a 
% given type, name and number of colors of the colorbrewer tables. 
% For more information on 'colorbrewer', please visit
% http://colorbrewer2.org/
%
% [colormap]=cbrewer(cname, ncol, interp_method)
%
% INPUT:
%  cname  name of colortable. One of the following
%      Sequential tables:
%        'Blues'
%        'BuGn'
%        'BuPu'
%        'GnBu'
%        'Greens'
%        'Greys'
%        'Oranges'
%        'OrRd'
%        'PuBu'
%        'PuBuGn'
%        'PuRd'
%        'Purples'
%        'RdPu'
%        'Reds'
%        'YlGn'
%        'YlGnBu'
%        'YlOrBr'
%        'YlOrRd'
%      Divergent tables:
%        'BrBG'
%        'PiYG'
%        'PRGn'
%        'PuOr'
%        'RdBu'
%        'RdGy'
%        'RdYlBu'
%        'RdYlGn'};
%      Qualitative tables:
%        'Accent'
%        'Dark2'
%        'Paired'
%        'Pastel1'
%        'Pastel2'
%        'Set1'
%        'Set2'
%        'Set3'
%   ncol  number of color in the table.
%   interp_method  if the table need to be interpolated, what method
%                  should be used for interp1? Default='pchip' (aka 'cubic')
% 
% A note on the number of colors: Based on the original data, there is
% only a certain number of colors available for each type and name of
% colortable. When 'ncol' is larger then the maximum number of colors
% originally given, an interpolation routine is called (interp1) to produce 
% the "extended" colormaps.
%
% See also: plot_brewer_cmap
%

% 
% The tables were generated from an MS-Excel file provided on the website
% http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html
%
% This product includes color specifications and designs developed by
% Cynthia Brewer (http://colorbrewer.org/).
% 
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 06.12.2011


% load colorbrewer data
load('colorbrewer.mat')
% initialise the colormap is there are any problems
colormap=[];
if (~exist('interp_method', 'var'))
    interp_method='pchip';
end


ctype_names={'div', 'seq', 'qual'};
ctype = [];
for try_ctype = ctype_names
  try_ctype = try_ctype{1};
  if isfield(colorbrewer.(try_ctype),cname)
    ctype = try_ctype;
    break;
  end
end
if isempty(ctype)
  error(['Could not find colormap named: ' cname]);
end

if (ncol>length(colorbrewer.(ctype).(cname)))
  cbrew_init=colorbrewer.(ctype).(cname){length(colorbrewer.(ctype).(cname))};
  colormap=interpolate_cbrewer(cbrew_init, interp_method, ncol);
  colormap=colormap./255;
  return
end

if isempty(colorbrewer.(ctype).(cname){ncol})
  
  while isempty(colorbrewer.(ctype).(cname){ncol})
    ncol=ncol+1;
  end        
  warning( ...
    ['The minimum number of colors for table ' cname ...
     ' is ' num2str(ncol) '.'])
end

colormap=(colorbrewer.(ctype).(cname){ncol})./255;

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
