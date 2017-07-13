function [varargout] = myaa(varargin)
%MYAA Render figure with anti-aliasing.
%   MYAA
%   Anti-aliased rendering of the current figure. This makes graphics look 
%   a lot better than in a standard matlab figure, which is useful for  
%   publishing results on the web or to better see the fine details in a
%   complex and cluttered plot. Some simple keyboard commands allow
%   the user to set the rendering quality interactively, zoom in/out and
%   re-render when needed.
%
%   Usage:
%     myaa: Renders an anti-aliased version of the current figure.
%
%     myaa(K): Sets the supersampling factor to K. Using a 
%     higher K yields better rendering but takes longer time. If K is 
%     omitted, it defaults to 4. It may be useful to run e.g. myaa(2) to 
%     make a low-quality rendering as a first try, because it is a lot 
%     faster than the default.
%
%     myaa([K D]): Sets supersampling factor to K but downsamples the 
%     image by a factor of D. If D is larger than K, the rendered image
%     will be smaller than the original. If D is smaller than K, the 
%     rendering will be bigger.
%
%     myaa('publish'): An experimental parameter, useful for publishing 
%     matlab programs (see example 3). Beware, it kills the current figure.
%
%   Interactivity:
%     The anti-aliased figure can be updated with the following keyboard 
%     commands:
%
%     <space>       Re-render image (to reflect changes in the figure)
%     +             Zoom in (decrease downsampling factor)
%     -             Zoom out (increase downsampling factor)
%     1 ... 9       Change supersampling and downsampling factor to ...
%     q             Quit, i.e. close the anti-aliased figure
%
%   Myaa can also be called with up to 3 parameters.
%   [FIG,RAW] = MYAA(K,AAMETHOD,FIGMODE)
%   [RAW] = MYAA('raw') Alec 2012
%
%   Parameters and output:
%     K         Subsampling factor. If a vector is specified, [K D], then 
%               the second element will describe the downsampling factor. 
%               Default is K = 4 and D = 4.
%     AAMETHOD  Downsampling method. Normally this is chosen automatically.
%               'standard': convolution based filtering and downsampling
%               'imresize': downsampling using the imresize command from
%               the image toolbox.
%               'noshrink': used internally
%     FIGMODE   Display mode
%               'figure': the normal mode, a new figure is created
%               'update': internally used for interactive sessions
%               'publish': used internally
%     FIG       A handle to the new anti-aliased figure
%
%     RAE       Raw image data (Alec 2012)
%     
%
%   Example 1:
%     spharm2;
%     myaa;
%     % Press '1', '2' or '4' and try '+' and '-'
%     % Press 'r' or <space> to update the anti-aliased rendering, e.g. 
%     % after rotating the 3-D object in the original figure.
%
%   Example 2:
%     line(randn(2500,2)',randn(2500,2)','color','black','linewidth',0.01)
%     myaa(8);
%
%   Example 3:
%     xpklein;
%     myaa(2,'standard');
%
%   Example 3:
%     Put the following in test.m
%        %% My test publish
%        %  Testing to publish some anti-aliased images
%        %
%        spharm2;          % Produce some nice graphics
%        myaa('publish');  % Render an anti-aliased version
%   Example 4:
%     spharm2;
%     imwrite(myaa('raw'),'myaa.png');
%
%     Then run:
%        publish test.m;
%        showdemo test;
%
%
%   BUGS:
%     Dotted and dashed lines in plots are not rendered correctly. This is
%     probably due to a bug in Matlab and it will hopefully be fixed in a
%     future version.
%     The OpenGL renderer does not always manage to render an image large
%     enough. Try the zbuffer renderer if you have problems or decrease the
%     K factor. You can set the current renderer to zbuffer by running e.g. 
%     set(gcf,'renderer','zbuffer').
%
%   See also PUBLISH, PRINT
%
%   Version 1.1, 2008-08-21
%   Version 1.0, 2008-08-05
%
%   Author: Anders Brun
%           anders@cb.uu.se
%

% This was adapted to produce an output image according to:
% http://www.alecjacobson.com/weblog/?p=2662#comment-9504
% 

%% Force drawing of graphics 
drawnow;

%% Find out about the current DPI...
screen_DPI = get(0,'ScreenPixelsPerInch');
self.scale = 1;

%% Determine the best choice of convolver.
% If IPPL is available, imfilter is much faster. Otherwise it does not 
% matter too much.
try
    if exists('imfilter','func')
        myconv = @imfilter;
    else
        myconv = @conv2;
    end
catch
    myconv = @conv2;
end

%% Set default options and interpret arguments
if isempty(varargin)
    self.K = [4 4];
    try
        imresize(zeros(2,2),1);
        self.aamethod = 'imresize';
    catch
        self.aamethod = 'standard';
    end
    self.figmode = 'figure';
elseif strcmp(varargin{1},'publish')
    self.K = [4 4];
    self.aamethod = 'noshrink';
    self.figmode = 'publish';
elseif strcmp(varargin{1},'update')
    self = get(gcf,'UserData');
    figure(self.source_fig);
    drawnow;
    self.figmode = 'update';
elseif strcmp(varargin{1},'lazyupdate')
    self = get(gcf,'UserData');
    self.figmode = 'lazyupdate';
elseif length(varargin) == 1
    if (ischar(varargin{1}) && strcmp(varargin{1},'raw')) || ...
      (iscell(varargin{1}) && strcmp(varargin{1}{1},'raw'))
      self.K = [4 4];
      self.figmode = 'raw';
      if iscell(varargin{1})
        self.scale = varargin{1}{2};
      end
    else
      self.K = varargin{1};
    self.figmode = 'figure';
    end
    if length(self.K) == 1
        self.K = [self.K self.K];
    end
    if self.K(1) > 16
        error('To avoid excessive use of memory, K has been limited to max 16. Change the code to fix this on your own risk.');
    end
    try
        imresize(zeros(2,2),1);
        self.aamethod = 'imresize';
    catch
        self.aamethod = 'standard';
    end
elseif length(varargin) == 2
    self.K = varargin{1};
    self.aamethod = varargin{2};
    self.figmode = 'figure';
elseif length(varargin) == 3
    self.K = varargin{1};
    self.aamethod = varargin{2};
    self.figmode = varargin{3};
    if strcmp(self.figmode,'publish') && ~strcmp(varargin{2},'noshrink')
        printf('\nThe AAMETHOD was not set to ''noshrink'': Fixed.\n\n');
        self.aamethod = 'noshrink';
    end
else
    error('Wrong syntax, run: help myaa');
end

if length(self.K) == 1
    self.K = [self.K self.K];
end

%% Capture current figure in high resolution
if ~strcmp(self.figmode,'lazyupdate');
    tempfile = 'myaa_temp_screendump.png';
    self.source_fig = gcf;
    current_paperpositionmode = get(self.source_fig,'PaperPositionMode');
    current_inverthardcopy = get(self.source_fig,'InvertHardcopy');
    set(self.source_fig,'PaperPositionMode','auto');
    set(self.source_fig,'InvertHardcopy','off');
    print(self.source_fig,['-r',num2str(self.scale*screen_DPI*self.K(1))], '-dpng', tempfile);
    set(self.source_fig,'InvertHardcopy',current_inverthardcopy);
    set(self.source_fig,'PaperPositionMode',current_paperpositionmode);
    self.raw_hires = imread(tempfile);
    delete(tempfile);

    %% Alec: This only works if the background color is white. It also seems to
    %% mess up grid axis lines etc.
    %self.raw_hires = print('-RGBImage',['-r',num2str(self.scale*screen_DPI*self.K(1))]);
end
%% Start filtering to remove aliasing
w = warning;
warning off;
if strcmp(self.aamethod,'standard') || strcmp(self.aamethod,'noshrink')
    % Subsample hires figure image with standard anti-aliasing using a
    % butterworth filter    
    kk = lpfilter(self.K(2)*3,self.K(2)*0.9,2);
    mm = myconv(ones(size(self.raw_hires(:,:,1))),kk,'same');
    a1 = max(min(myconv(single(self.raw_hires(:,:,1))/(255),kk,'same'),1),0)./mm;
    a2 = max(min(myconv(single(self.raw_hires(:,:,2))/(255),kk,'same'),1),0)./mm;
    a3 = max(min(myconv(single(self.raw_hires(:,:,3))/(255),kk,'same'),1),0)./mm;
    if strcmp(self.aamethod,'standard')
        if abs(1-self.K(2)) > 0.001
            raw_lowres = double(cat(3,a1(2:self.K(2):end,2:self.K(2):end),a2(2:self.K(2):end,2:self.K(2):end),a3(2:self.K(2):end,2:self.K(2):end)));
        else
            raw_lowres = self.raw_hires;
        end
    else
        raw_lowres = double(cat(3,a1,a2,a3));
    end
elseif strcmp(self.aamethod,'imresize')
    % This is probably the fastest method available at this moment...
    raw_lowres = single(imresize(self.raw_hires,1/self.K(2),'bilinear'))/255;
end
warning(w);

%% Place the anti-aliased image in some image on the screen ...
if strcmp(self.figmode,'figure');
    % Create a new figure at the same place as the previous
    % The content of this new image is just a bitmap...
    oldpos = get(gcf,'Position');
    self.myaa_figure = figure;
    fig = self.myaa_figure;
    set(fig,'Menubar','none');
    set(fig,'Resize','off');
    sz = size(raw_lowres);
    set(fig,'Units','pixels');
    pos = [oldpos(1:2) sz(2:-1:1)];
    set(fig,'Position',pos);
    ax = axes;
    image(raw_lowres);
    set(ax,'Units','pixels');
    set(ax,'Position',[1 1 sz(2) sz(1)]);
    axis off;
elseif strcmp(self.figmode,'publish');
    % Create a new figure at the same place as the previous
    % The content of this new image is just a bitmap...
    self.myaa_figure = figure;
    fig = self.myaa_figure;
    current_units = get(self.source_fig,'Units');
    set(self.source_fig,'Units','pixels');
    pos = get(self.source_fig,'Position');
    set(self.source_fig,'Units',current_units);
    set(fig,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    ax = axes;
    image(raw_lowres);
    set(ax,'Units','normalized');
    set(ax,'Position',[0 0 1 1]);
    axis off;
    close(self.source_fig);
elseif strcmp(self.figmode,'update');
    fig = self.myaa_figure;
    figure(fig);
    clf;    
    set(fig,'Menubar','none');
    set(fig,'Resize','off');
    sz = size(raw_lowres);
    set(fig,'Units','pixels');
    pos = get(fig,'Position');
    pos(3:4) = sz(2:-1:1);
    set(fig,'Position',pos);
    ax = axes;
    image(raw_lowres);
    set(ax,'Units','pixels');
    set(ax,'Position',[1 1 sz(2) sz(1)]);
    axis off;
elseif strcmp(self.figmode,'lazyupdate');
    clf;
    fig = self.myaa_figure;
    sz = size(raw_lowres);
    pos = get(fig,'Position');
    pos(3:4) = sz(2:-1:1);
    set(fig,'Position',pos);
    ax = axes;
    image(raw_lowres);
    set(ax,'Units','pixels');
    set(ax,'Position',[1 1 sz(2) sz(1)]);
    axis off;    
elseif strcmp(self.figmode,'raw')
  % don't create a new figure
end

%% Store current state

set(gcf,'userdata',self);
set(gcf,'KeyPressFcn',@keypress);
set(gcf,'Interruptible','off');

%% Avoid unnecessary console output
if strcmp(self.figmode,'raw')
  varargout{1} = raw_lowres;
else
  if nargout == 1
    varargout{1} = fig;
  elseif nargout ==2
    varargout{1} = fig;
    varargout{2} = raw_lowres;
  end
end

%% A simple lowpass filter kernel (Butterworth).
% sz is the size of the filter
% subsmp is the downsampling factor to be used later
% n is the degree of the butterworth filter
function kk = lpfilter(sz, subsmp, n)
sz = 2*floor(sz/2)+1; % make sure the size of the filter is odd
cut_frequency = 0.5 / subsmp;
range = (-(sz-1)/2:(sz-1)/2)/(sz-1);
[ii,jj] = ndgrid(range,range);
rr = sqrt(ii.^2+jj.^2);
kk = ifftshift(1./(1+(rr./cut_frequency).^(2*n)));
kk = fftshift(real(ifft2(kk)));
kk = kk./sum(kk(:));

function keypress(src,evnt)
if isempty(evnt.Character)
    return
end
recognized = 0;
self = get(gcf,'userdata');

if evnt.Character == '+'
    self.K(2) = max(self.K(2).*0.5^(1/2),1);    
    recognized = 1;
    set(gcf,'userdata',self);
    myaa('lazyupdate');
elseif evnt.Character == '-'
    self.K(2) = min(self.K(2).*2^(1/2),16);
    recognized = 1;
    set(gcf,'userdata',self);
    myaa('lazyupdate');
elseif evnt.Character == ' ' || evnt.Character == 'r' || evnt.Character == 'R'
    set(gcf,'userdata',self);
    myaa('update');
elseif evnt.Character == 'q' 
    close(gcf); 
elseif find('123456789' == evnt.Character)
    self.K = [str2double(evnt.Character) str2double(evnt.Character)];
    set(gcf,'userdata',self);
    myaa('update');
end


% Copyright (c) 2009, Anders Brun
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
