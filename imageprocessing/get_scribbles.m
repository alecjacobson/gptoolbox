function [SI,S,C] = get_scribbles(im,varargin)
  % Display an image and ask the user to draw scribbles then hit enter
  %
  % [SI,S,C] = get_scribbles(im)
  %
  % Inputs:
  %   im  w by h by channels image
  %   Optional:
  %     'LineWidth' followed by line width of brush
  % Outputs:
  %   SI  w by h scribbles mask (0 no scribble, otherwise scribble id)
  %   S  w by h scribbles mask in RGB
  %   C  colors corresponding to each scribble
  %

  f = gcf;
  h = imshow(im);

  line_width = 10;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'LineWidth'},{'line_width'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  rotate3d off;
  set(h,'ButtonDownFcn',@ondown);
  set(gcf,'keypressfcn',        @onkeypress);
  t = title('Draw scribbles, then hit ENTER');

  P = {};
  p = {};
  pc = 0;

  done = false;
  while(~done)
    drawnow;
  end
  
  %set(h,'CData',0*get(h,'CData'));
  %f_pos = get(gcf,'Position');
  %f_pos(3:4) = [size(im,2) size(im,1)];
  %set(gcf,'Position',f_pos);
  %set(gcf,'Color', [0,0,0]);
  %set(gca,'visible', 'off');
  %set(gca,'Position',[0 0 1 1]);
  old_color = get(gcf,'Color');
  set(h,'CData',0*get(h,'CData'));
  set(t,'Visible','off');
  set(gcf,'Color',[1 0 1]);
  old_aa = get(gcf,'GraphicsSmoothing');
  set(gcf,'GraphicsSmoothing','off')
  drawnow;
  F = getframe(gcf);
  S = F.cdata;
  S = imresize(imtrim(im2double(S)),[size(im,1) size(im,2)],'nearest');
  set(gcf,'Color',old_color);
  set(gcf,'GraphicsSmoothing',old_aa);
  %S = S(1:size(im,1),1:size(im,2),:);
  SI = rgb2ind(S,pc+1);
  if(numel(unique(SI(:))) ~= pc+1)
    warning(['IDs not right']);
  end

  % match class of input
  switch class(im)
  case 'double'
    S = im2double(S);
  case 'single'
    S = im2single(S);
  case 'uint8'
    S = im2uint8(S);
  case 'uint16'
    S = im2uint16(S);
  case 'int16'
    S = im2int16(S);
  case 'int16'
    S = im2int16(S);
  otherwise
    warning(['Input is strange class: ' class(im)]);
  end

  % Callback for mouse press
  function ondown(src,ev)
    % Tell window that we'll handle drag and up events
    set(gcf,'windowbuttonmotionfcn', @ondrag);
    set(gcf,'windowbuttonupfcn',     @onup);
    % increment plot count
    pc=pc+1;
    append_current_point();
  end

  % Callback for mouse press
  function ondrag(src,ev)
    append_current_point();
  end

  % Callback for mouse release
  function onup(src,ev)
    % Tell window to handle drag and up events itself
    set(gcf,'windowbuttonmotionfcn','');
    set(gcf,'windowbuttonupfcn','');
  end

  function onkeypress(src,ev)
    % escape character id
    ESC = char(27);
    ENTER = char(13);
    switch ev.Character
    case ESC
      finish();
    case ENTER
      finish();
    otherwise
      warning(['Unknown key: ' ev.Character ...
        ' (' num2str(uint8(ev.Character(1))) ')']);
    end
  end

  function append_current_point()
    % get current mouse position, and remember old one
    cp = get(gca,'currentpoint');
    if pc > numel(P)
      P{pc} = [];
    end
    P{pc} = [P{pc};cp(1,:)];
    if pc > numel(p) || isempty(p{pc})
      hold on;
      p{pc} = plot(P{pc}(:,1),P{pc}(:,2),'LineWidth',line_width,'Color',next_color);
      hold off;
    else
      set(p{pc},'Xdata',P{pc}(:,1),'Ydata',P{pc}(:,2));
    end
  end

  function nc = next_color()
    if (pc-1) <= 20
      nc = [(pc-1)/20 0 1];
    else
      nc = [1 (pc-1)/40 1];
    end
  end

  function finish()
    done = true;
    set(gcf,'windowbuttonmotionfcn','');
    set(gcf,'windowbuttonupfcn','');
    set(h,'ButtonDownFcn','');
    set(gcf,'keypressfcn','');
  end
end
