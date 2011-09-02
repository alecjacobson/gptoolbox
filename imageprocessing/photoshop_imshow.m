function photoshop_imshow(im,z)
  % PHOTOSHOP_IMSHOW show an image like photoshop with a given zoom. If zoomed in
  % enough then faint lines appear separating pixels
  %
  % photoshop_imshow(im,z)
  %
  % Inputs:
  %   im  image to show
  %   z  zoom amount, 
  %
  % Note sure how this will work if the plot is zoomed, rotated or put in a
  % subplot


  % current figure handle
  fh = gcf;

  if ~exist('z','var')
    z = 1;
  end

  % zoom scale should be positive
  assert(z > 0);
  % zoom scale should be either smaller than one or an integer
  if z > 1
    z = round(z);
  end

  % remove any resize function if one already exists
  set(fh,'ResizeFcn','');

  % resize image
  rim = photoshop_imresize(im,z);
  ih = imshow(rim);
  % get handle to axes
  ah = gca;

  set(fh,'ResizeFcn',@onresize);
  zh = zoom(fh);
  % keep track of range of x
  set(zh,'ButtonDownFilter',@zoomfilter);
  %set(zh,'ActionPreCallback',@beforezoom);
  %set(zh,'ActionPostCallback',@afterzoom);

  AP = get(ah,'Position');
  %xlim = get(ah,'XLim');
  %ylim = get(ah,'YLim');

  %P = get(gcf,'Position');
  %set(gcf,'Position',[P(1:2) z*size(im,2) z*size(im,1)]);
  %set(gca,'Position',[0 0 1 1]);

  function onresize(src,ev)
    %% axis position
    %AP = get(ah,'Position');
    % figure position
    P = get(fh,'Position');
    % displayed scale of image is min of axis position times figure position
    %d = min((AP([4 3]).*P([4 3]))./[size(im,1) size(im,2)]);
    z = min((AP([4 3]).*P([4 3]))./[size(im,1) size(im,2)]);
    if z > 1
      z = round(z);
    end
    z = max([z 1/min([size(im,1) size(im,2)])]);
    rim = photoshop_imresize(im,z);
    assert(min(size(rim)) > 0)
    % temporarily turn off resize func to prevent recursion
    set(fh,'ResizeFcn','');
    set(ih,'CData',rim);
    set(ih,'XData',[1 size(rim,2)]);
    set(ih,'YData',[1 size(rim,1)]);
    set(ah,'XLim',[0.5 size(rim,2)+0.5]);
    set(ah,'YLim',[0.5 size(rim,1)+0.5]);
    %xlim = get(ah,'XLim');
    %ylim = get(ah,'YLim');

    WH = [size(rim,2) size(rim,1)]./P([3 4]);
    set(ah,'Position',[(1-WH)*0.5 WH]);
    set(fh,'ResizeFcn',@onresize);
  end


  % zoom filter
  function res = zoomfilter(src,ev)
    res = false;
  end

  %prev_xlim = [];
  %prev_ylim = [];
  %% pre zoom callback
  %function beforezoom(src,ev)
  %  prev_xlim = get(ah,'XLim');
  %  prev_ylim = get(ah,'YLim');
  %  get(ah,'Xlim')
  %end

  %zz = 1;
  %% post zoom callback
  %function afterzoom(src,ev)
  %  get(ah,'Xlim')
  %  new_xlim = get(ah,'XLim');
  %  new_ylim = get(ah,'YLim');
  %  % zoom scale
  %  zz = zz*sum(prev_xlim.*[-1 1])/sum(new_xlim.*[-1 1]);
  %  % temporarily turn off resize func to prevent recursion
  %  set(fh,'ResizeFcn','');
  %  rim = photoshop_imresize(im,zz*z);
  %  set(ih,'CData',rim);
  %  %set(ah,'Xlim',(new_xlim-xlim(1))*zz+new_xlim(1));
  %  %set(ah,'Ylim',(new_ylim-ylim(1))*zz+new_ylim(1));

  %  set(ih,'XData',[1 size(rim,2)]);
  %  set(ih,'YData',[1 size(rim,1)]);
  %  set(ah,'Xlim',(new_xlim-0.5)*zz+0.5)
  %  set(ah,'Ylim',(new_ylim-0.5)*zz+0.5);
  %  %set(ah,'XLim',[0.5 size(rim,2)+0.5]);
  %  %set(ah,'YLim',[0.5 size(rim,1)+0.5]);
  %  %xlim = get(ah,'XLim');
  %  %ylim = get(ah,'YLim');
  %  %%%imshow(rim,'Parent',ah);
  %  %rxlim = get(ah,'Xlim');
  %  %rylim = get(ah,'Ylim');
  %  %set(ah,'Xlim',(new_xlim-xlim(1))*zz+rxlim(1));
  %  %set(ah,'Ylim',(new_ylim-ylim(1))*zz+rylim(1));
  %  set(fh,'ResizeFcn',@onresize);
  %end

end
