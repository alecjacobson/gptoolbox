function tilefigs(tile,border)
% <cpp> tile figure windows usage: tilefigs ([nrows ncols],border_in pixels)
% nrows,ncols may be nonpositive or inf to mean that they should be computed
% automatically
% Restriction: maximum of 100 figure windows
% Without arguments, tilefigs will determine the closest N x N grid


%Charles Plum                    Nichols Research Corp.
%<cplum@nichols.com>             70 Westview Street
%Tel: (781) 862-9400             Kilnbrook IV
%Fax: (781) 862-9485             Lexington, MA 02173


maxpos  = get (0,'screensize'); % determine terminal size in pixels
maxpos(4) = maxpos(4) - 25;
hands   = get (0,'Children');   % locate fall open figure handles
hands   = sort(hands);          % sort figure handles
numfigs = size(hands,1);        % number of open figures


maxfigs = 100;


if (numfigs>maxfigs)            % figure limit check
  disp([' More than ' num2str(maxfigs) ' figures ... get serious pal'])
  return
end


if nargin > 0 
  nrows = tile(1);
  ncols = tile(2);
  both_inf = false;
  
  if (ncols <= 0 || ncols == inf) && (nrows <= 0 || nrows == inf)
    both_inf = true;
  elseif ncols <= 0 || ncols == inf
      ncols = ceil(numfigs/nrows);
  elseif nrows <= 0 || nrows == inf
      nrows = ceil(numfigs/ncols);
  end
  
  if numfigs > nrows*ncols
    disp ([' requested tile size too small for ' ...
        num2str(numfigs) ' open figures '])
        return
  end
end

if (nargin == 0 || both_inf)
  maxfactor = sqrt(maxfigs);       % max number of figures per row or column
  sq = [2:maxfactor].^2;           % vector of integer squares
  sq = sq(find(sq>=numfigs));      % determine square grid size
  gridsize = sq(1);                % best grid size
  nrows = sqrt(gridsize);          % figure size screen scale factor
  ncols = nrows;                   % figure size screen scale factor
end

if nargin < 2
  border = 0;
else
  maxpos(3) = maxpos(3) - 2*border;
  maxpos(4) = maxpos(4) - 2*border;
end
xlen = fix(maxpos(3)/ncols) - 30; % new tiled figure width
ylen = fix(maxpos(4)/nrows) - 45; % new tiled figure height


% tile figures by postiion 
% Location (1,1) is at bottom left corner
pnum=0;
for iy = 1:nrows
  ypos = maxpos(4) - fix((iy)*maxpos(4)/nrows) + border +25; % figure location (row)
  for ix = 1:ncols
        xpos = fix((ix-1)*maxpos(3)/ncols + 1) + border+7;     % figure location (column)
        pnum = pnum+1;
    if (pnum>numfigs)
                break
        else
          figure(hands(pnum))
      set(hands(pnum),'Position',[ xpos ypos xlen ylen ]); % move figure
        end
  end
end
return


% -------------------------------------------------------------------------
