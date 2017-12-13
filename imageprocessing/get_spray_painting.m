function [im,sh] = get_spray_painting()
  % GET_SPRAY_PAINTING Interactively paint on top of the current image
  % 
  % Outputs:
  %   im  image
  %
  figure(gcf);
  set(gcf,'windowbuttondownfcn',    @ondown);
  set(gcf,'KeyPressFcn',          @onkeypress    );
  hold on;
  A = axis();
  % Initialize with grayscale nan image
  h = round(A(4)-A(1));
  w = round(A(2)-A(1));
  ratio = h/w;
  if h <= 1 || w <=1 
    w = 512;
    h = ceil(ratio*w);
  end
  im = ones(h,w);
  [X,Y] = meshgrid([A(1) A(2)],[A(3) A(4)]);
  sh = surf(X,Y,1+0*X,'CData',im,'FaceColor','texturemap');
  hold off;
  drawnow;
  last_pos = [];
  pos = [];
  % window/brush radius
  r = 20;
  K = matrixnormalize(fspecial('gaussian',r*2+1,5));

  still_painting = true;
  while still_painting
    drawnow;
  end

  function paint_at(p)
    p = round( ((p-[X(1) Y(1)])./([X(end) Y(end)]-[X(1) Y(1)])).*[w h]);
    I = max(1,p(2)-r):min(size(im,1),p(2)+r);
    J = max(1,p(1)-r):min(size(im,2),p(1)+r);
    Kij = K(I-p(2)+r+1,J-p(1)+r+1);
    im(I,J) = min(max(im(I,J)-Kij,0),1);
    sh.CData = im;
  end
  function ondown(src,ev)
    set(gcf,'windowbuttonmotionfcn',@ondrag);
    set(gcf,'windowbuttonupfcn',    @onup);
    last_pos = eye(1,2)*get(gca,'currentpoint')*eye(3,2);
    pos = last_pos;
    paint_at(pos);
  end
  function ondrag(src,ev)
    %% paint from last_pos to pos
    last_pos = pos;
    pos = eye(1,2)*get(gca,'currentpoint')*eye(3,2);
    paint_at(pos);
  end
  function onup(src,ev)
    set(gcf,'windowbuttonmotionfcn','');
    set(gcf,'windowbuttonupfcn','');
    last_pos = [];
    pos = [];
  end
  function onkeypress(src,ev)
    %switch ev.Character
    %case 'r'
    %end
    still_painting = false;
  end
end
