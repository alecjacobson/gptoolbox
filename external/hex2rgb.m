function [ rgb ] = hex2rgb(hex,range)
% hex2rgb converts hex color values to rgb arrays on the range 0 to 1. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% SYNTAX:
% rgb = hex2rgb(hex) returns rgb color values in an n x 3 array. Values are
%                    scaled from 0 to 1 by default. 
%                    
% rgb = hex2rgb(hex,256) returns RGB values scaled from 0 to 255. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% EXAMPLES: 
% 
% myrgbvalue = hex2rgb('#334D66')
%    = 0.2000    0.3020    0.4000
% 
% 
% myrgbvalue = hex2rgb('334D66')  % <-the # sign is optional 
%    = 0.2000    0.3020    0.4000
% 
%
% myRGBvalue = hex2rgb('#334D66',256)
%    = 51    77   102
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myrgbvalues = hex2rgb(myhexvalues)
%    =   0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%        0.8000    0.6000    0.2000
%        0.2000    0.2000    0.9020
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myRGBvalues = hex2rgb(myhexvalues,256)
%    =   51    77   102
%       128   153   179
%       204   153    51
%        51    51   230
% 
% HexValsAsACharacterArray = {'#334D66';'#8099B3';'#CC9933';'#3333E6'}; 
% rgbvals = hex2rgb(HexValsAsACharacterArray)
% 
% * * * * * * * * * * * * * * * * * * * * 
% Chad A. Greene, April 2014
%
% Updated August 2014: Functionality remains exactly the same, but it's a
% little more efficient and more robust. Thanks to Stephen Cobeldick for
% the improvement tips. In this update, the documentation now shows that
% the range may be set to 256. This is more intuitive than the previous
% style, which scaled values from 0 to 255 with range set to 255.  Now you
% can enter 256 or 255 for the range, and the answer will be the same--rgb
% values scaled from 0 to 255. Function now also accepts character arrays
% as input. 
% 
% * * * * * * * * * * * * * * * * * * * * 
% See also rgb2hex, dec2hex, hex2num, and ColorSpec. 
% 

%% Input checks:

assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.') 

if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end

%% Tweak inputs if necessary: 

if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
    
    % In case cell array elements are separated by a comma instead of a
    % semicolon, reshape hex:
    if isrow(hex)
        hex = hex'; 
    end
    
    % If input is cell, convert to matrix: 
    hex = cell2mat(hex);
end

if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end

if nargin == 1
    range = 1; 
end

%% Convert from hex to rgb: 

switch range
    case 1
        rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;

    case {255,256}
        rgb = reshape(sscanf(hex.','%2x'),3,[]).';
    
    otherwise
        error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end

end

