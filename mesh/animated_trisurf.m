function h = animated_trisurf(F,X,Y,Z)
  % ANIMATED_TRISURF  animated, interactive trisurf plot
  %
  % h = animated_trisurf(F,X,Y,Z)
  %
  % Inputs:
  %   F  #F by 3 by {#frames | 1} list of triangle indices
  %   V  #V by dim by #{frames | 1} list of vertex positions
  %    or
  %   X  #V by #{frames | 1} list of X positions
  %   Y  #V by #{frames | 1} list of Y positions
  %   Z  #V by #{frames | 1} list of Z positions
  % Outputs:
  %   h  handle to plot
  %
  % See also: tsurf, trisurf
  %

  if nargin == 2
    V = X;
    dim = size(V,2);
    X = squeeze(V(:,1,:));
    Y = squeeze(V(:,2,:));
    if dim == 3
      Z = squeeze(V(:,3,:));
    else
      Z = 0*X(:,1);
    end
  end

  
  % number of mesh vertices
  n = size(X,1);
  assert(size(Y,1) == n);
  assert(size(Z,1) == n);

  % count of frames
  fc = max([size(X,2) size(Y,2) size(Z,2)]);
  assert(size(X,2) == 1 || size(X,2) == fc);
  assert(size(Y,2) == 1 || size(Y,2) == fc);
  assert(size(Z,2) == 1 || size(Z,2) == fc);

  assert(size(F,3) == 1 || size(F,3) == fc);

  % set animating to false
  animating = false;
  % current frame
  fi = 1;
  % current frame delta
  fd = 1;
  h = trisurf(F(:,:,fi),X(:,fi),Y(:,fi),Z(:,fi), ...
    'FaceColor','interp', ...
    'ButtonDownFcn',@onmeshdown);
  was_hold = ishold;
  hold on;
  p = plot3([max(X(:));min(X(:))],[max(Y(:));min(Y(:))],[max(Z(:));min(Z(:))]);
  if ~was_hold
    hold off;
  end
  set(p,'Visible','off');

  t = timer('TimerFcn',@increment_if_animating,'Period',0.03,'ExecutionMode','fixedDelay');
  start(t);


  % handle simple case where there's only a single frame
  if fc == 1
    return
  end


  % on mouse down on mesh
  function onmeshdown(src,ev)
    if(strcmp('normal',get(gcf,'SelectionType')))
      % left-click
      down_type = 'left';
    else
      % other (right) click
      down_type = 'right';
    end

    switch down_type
    case 'left'
      animating = ~animating;
    end
  end

  function increment_if_animating(src,ev)
    if animating
      increment();
    elseif ~ishandle(h)
      stop(t);
      delete(t);
      return;
    end
  end

  function increment()
    % if plot is closed then delete timer
    if ~ishandle(h)
      stop(t);
      delete(t);
      return;
    end

    fi = fi+fd;
    if fi <= 1
      fd = 1;
      fi = 1;
      animating = false;
    elseif fi >= fc
      fd = -1;
      fi = fc;
      animating = false;
    end

    % always check that handle still exists before updating
    if ishandle(h)
      set(h,'Vertices',[ ...
        X(:,mod(fi-1,size(X,2))+1) ...
        Y(:,mod(fi-1,size(Y,2))+1) ...
        Z(:,mod(fi-1,size(Z,2))+1)]);
      if size(F,3) == fc
        set(h,'Faces',F(:,:,mod(fi-1,fc)+1));
      end
      set(h,'CData',Z(:,mod(fi-1,size(Z,2))+1));
      drawnow
    end
  end


end
