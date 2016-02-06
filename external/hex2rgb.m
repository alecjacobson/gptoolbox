function RGB = hex2rgb(HEX)
% HEX2RGB - convert hexadecimal color strings to RGB values
%
%   RGB = HEX2RGB(HEX) converts the hexadimal color string HEX to its
%   corresponding RGB values. RGB has three columns representing the red,
%   green and blue component of the color.  
%   For a cell array of color strings, RGB will have as many rows as
%   elements of the cell array. For a character array HEX, RGB will have as
%   many rows as HEX. 
%
%   Three-digit hexadecimal color strings are expanded to six-digit strings
%   by doubling each digit (i.e., XYZ -> XXYYZZ).
%   
%   Examples:
%     hex2rgb('556b2f') % 6 digit string
%       % -> [ 85 107  47]
%     hex2rgb('f0f')    % 3 digit string  
%       % -> [255   0 255] 
%     hex2rgb({'8B4513','FF0'})  % cell array
%       % -> [139 69 19 ; 255 255 0]
%     hex2rgb(['FF6347' ; '40E0D0']) % character array with multiple rows
%       % -> [255 99 71 ; 64 224 208]
%
%   Hexadecimal color strings are three-byte triplets representing the red,
%   green and blue component. One byte represents a number in the range 00
%   to FF (in hexadecimal notation). More information:
%   http://en.wikipedia.org/wiki/Web_colors
%
%   See also HEX2DEC

% Created in Matlab R2011b 
% version 1.0 (feb 2014)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% 1.0 (feb 2014) - created

% check the input
if ~iscellstr(HEX)
    if ischar(HEX)
       HEX = cellstr(HEX) ;
    else
       error('Input should be a cell array of strings, or a character matrix.') ; 
    end
end

HEX = upper(HEX) ;

tf = cellfun(@(S) numel(S)==3, HEX) ;
if any(~tf) && ~all(cellfun(@(S) numel(S)==6, HEX(~tf)))
    error('All hexadecimal color strings should be 3 or 6 characters long.') ;
elseif any(tf)
    % convert any 3 character arrays to 6 character arrays
    % xyz -> xxyyzz
    HEX(tf) = cellfun(@(S) S([1 1 2 2 3 3]), HEX(tf),'un',0) ;
end

% convert to numbers between 0 and 15
HEX = double(char(HEX)) ;
tf = HEX > 64 ; % letters
HEX(~tf) = HEX(~tf) - '0' ;     % numbers:  0-9
HEX(tf)  = HEX(tf) - 'A' + 10 ; % letters: 10-15

if any(HEX(:) < 0) || any(HEX(:) > 15)
    error('Input string found with characters other than 0-9, a-f, or A-F.') ;
end

% reshape in two rows, each column is a byte (xy)
% the first 3 columns are the RGB values of the first column, etc.
HEX = reshape(HEX.',2,[]) ;
% [x ; y] -> [16*x ; y] 
HEX(1,:) = 16 * HEX(1,:) ;
% -> 16*x+y
RGB = sum(HEX,1) ; 

% reshape in RGB triplet
RGB = reshape(RGB,3,[]) ; % color triplets
RGB = RGB.' ; % we want color as rows




% Copyright (c) 2014, Jos (10584)
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
%     * Neither the name of the  nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
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
