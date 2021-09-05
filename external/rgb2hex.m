function [ hex ] = rgb2hex(rgb)
% rgb2hex converts rgb color values to hex color format. 
% 
% This function assumes rgb values are in [r g b] format on the 0 to 1
% scale.  If, however, any value r, g, or b exceed 1, the function assumes
% [r g b] are scaled between 0 and 255. 
% 
% * * * * * * * * * * * * * * * * * * * * 
% SYNTAX:
% hex = rgb2hex(rgb) returns the hexadecimal color value of the n x 3 rgb
%                    values. rgb can be an array. 
% 
% * * * * * * * * * * * * * * * * * * * * 
% EXAMPLES: 
% 
% myhexvalue = rgb2hex([0 1 0])
%    = #00FF00
% 
% myhexvalue = rgb2hex([0 255 0])
%    = #00FF00
% 
% myrgbvalues = [.2 .3 .4;
%                .5 .6 .7; 
%                .8 .6 .2;
%                .2 .2 .9];
% myhexvalues = rgb2hex(myrgbvalues) 
%    = #334D66
%      #8099B3
%      #CC9933
%      #3333E6
% 
% * * * * * * * * * * * * * * * * * * * * 
% Chad A. Greene, April 2014
% 
% Updated August 2014: Functionality remains exactly the same, but it's a
% little more efficient and more robust. Thanks to Stephen Cobeldick for
% his suggestions. 
% 
% * * * * * * * * * * * * * * * * * * * * 
% See also hex2rgb, dec2hex, hex2num, and ColorSpec. 

%% Check inputs: 

assert(nargin==1,'This function requires an RGB input.') 
assert(isnumeric(rgb)==1,'Function input must be numeric.') 

sizergb = size(rgb); 
assert(sizergb(2)==3,'rgb value must have three components in the form [r g b].')
assert(max(rgb(:))<=255& min(rgb(:))>=0,'rgb values must be on a scale of 0 to 1 or 0 to 255')

%% If no value in RGB exceeds unity, scale from 0 to 255: 
if max(rgb(:))<=1
    rgb = round(rgb*255); 
else
    rgb = round(rgb); 
end

%% Convert (Thanks to Stephen Cobeldick for this clever, efficient solution):

hex(:,2:7) = reshape(sprintf('%02X',rgb.'),6,[]).'; 
hex(:,1) = '#';


end
