function [SI,S,C] = get_scribbles(im)
  % Display an image and ask the user to draw scribbles then hit enter
  %
  % [SI,S,C] = get_scribbles(im)
  %
  % Inputs:
  %   im  w by h by channels image
  % Outputs:
  %   SI  w by h scribbles mask (0 no scribble, otherwise scribble id)
  %   S  w by h scribbles mask in RGB
  %   C  colors corresponding to each scribble
  %

  f = gcf;
  h = imshow(im);

  rotate3d off;
  set(h,'ButtonDownFcn',@ondown);
  set(gcf,'keypressfcn',        @onkeypress);
  title('Draw scribbles, then hit ENTER');

  P = {};
  p = {};
  pc = 0;

  done = false;
  while(~done)
    drawnow;
  end
  
  set(h,'CData',0*get(h,'CData'));
  drawnow;
  waitfor(h,'CData',0*get(h,'CData'));
  f_pos = get(gcf,'Position');
  f_pos(3:4) = [size(im,2) size(im,1)];
  set(gcf,'Position',f_pos);
  drawnow;
  set(gcf, 'Color', [0,0,0]);
  drawnow;
  waitfor(gcf, 'Color', [0,0,0]);
  set(gca, 'visible', 'off');
  drawnow;
  waitfor(gca, 'visible', 'off');
  set(gca,'Position',[0 0 1 1]);
  drawnow;
  waitfor(gca,'Position',[0 0 1 1]);
  drawnow;
  F = getframe(gca);
  drawnow;
  F = getframe(gca);
  drawnow;
  F = getframe(gca);
  S = F.cdata;
  S = S(1:size(im,1),1:size(im,2),:);
  [~,~,SI] = unique([0;reshape(rgb2gray(S),[],1)]);
  SI = reshape(SI(2:end)-1,size(S,1),size(S,2));
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
      p{pc} = plot(P{pc}(:,1),P{pc}(:,2),'LineWidth',2,'Color',next_color);
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
