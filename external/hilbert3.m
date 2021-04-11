function [x,y,z] = hilbert3(n)
% Hilbert 3D curve.
%
% function [x,y,z] = hilbert3(n) gives the vector coordinates of points
% in n-th order Hilbert curve of area 1.
%
% Example: plot the 3-rd order curve
%
% [x,y,z] = hilbert3(3); plot3(x,y,z)
%
% Inputs:
%   n  number of subdivision iterations 
% Outputs:
%   x  list of x coordinates
%   y  list of y coordinates
%   z  list of z coordinates
%
% Example:
%   n=3;
%   [x,y,z] = hilbert3(n);
%   P = [x;y;z]';
%   [V,F] = wire_mesh(P,[1:size(P,1)-1;2:size(P,1)]','Thickness',2^(1-n)/6);
%   [V,F] = loop(V,F,2);
%   tsurf(F,V);
%



%   Copyright (c) by Ivan Martynov
%   Inspired by function [x,y] = hilbert(n) made by Federico Forte
%   Date: September 17, 2009

if nargin ~= 1
    n = 2;
end

if n <= 0
    x = 0;
    y = 0;
    z = 0;
else
    [xo,yo,zo] = hilbert3(n-1);
    x = .5*[.5+zo .5+yo -.5+yo -.5-xo -.5-xo -.5-yo .5-yo .5+zo];
    y = .5*[.5+xo .5+zo .5+zo .5+yo -.5+yo -.5-zo -.5-zo -.5-xo];
    z = .5*[.5+yo -.5+xo -.5+xo .5-zo .5-zo -.5+xo -.5+xo .5-yo];
end

%{
Copyright (c) 2009, Ivan Martynov
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
%}
