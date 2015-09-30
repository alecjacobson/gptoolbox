function f = plot_brewer_cmap()
% Plots and identifies the various colorbrewer tables available.
% Is called by cbrewer.m when no arguments are given.
% 
% f = plot_brewer_cmap()
%
% Outputs:
%  f  handle to new figure
%
% See also: cbrewer
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 14.10.2011
%

  load('colorbrewer.mat')
  
  ctypes={'div', 'seq', 'qual'};
  ctypes_title={'div', 'seq', 'qual'};
  cnames{1,:}={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};
  cnames{2,:}={'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
               'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd'};
  cnames{3,:}={'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'};
  
  f = figure('position', [314 327 807 420])
  for itype=1:3
      
      %fh(itype)=figure();
      
      subplot(1,3,itype)
      
      for iname=1:length(cnames{itype,:})
          
          ncol=length(colorbrewer.(ctypes{itype}).(cnames{itype}{iname}));
          fg=1./ncol; % geometrical factor
  
          X=fg.*[0 0 1 1];
          Y=0.1.*[1 0 0 1]+(2*iname-1)*0.1;
          F=cbrewer(cnames{itype}{iname}, ncol);
  
          for icol=1:ncol
              X2=X+fg.*(icol-1);
              fill(X2,Y,F(icol, :), 'linestyle', 'none')
              text(-0.1, mean(Y), cnames{itype}{iname}, 'HorizontalAlignment', 'right', 'FontSize', 18);
              xlim([-0.4, 1])
              hold all
          end % icol
          %set(gca, 'box', 'off')
          title(ctypes_title{itype}, 'FontWeight', 'bold', 'FontSize', 16, 'FontName' , 'AvantGarde')
          axis off
          set(gcf, 'color', [1 1 1])
      end % iname
  
  end %itype
  
  set(gcf, 'MenuBar', 'none')
  set(gcf, 'Name', 'ColorBrewer Color maps')
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
