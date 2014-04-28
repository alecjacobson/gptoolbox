function deform_mls(varargin)
  % DEFORM
  %
  % deform(V,F,C,mlsd)
  %
  % Deform a mesh (V,F) using linear blend skinning, by applying
  % transformations at a set of control points C and propogating the
  % deformation to the mesh via MLS data mlsd
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices
  %  C  list of control point positions
  %  mlsd  MLS data struct
  %  Options:
  %    'Lines' followed by a #L by 2 list of indices into C for the line
  %      handles. Means that mlsd.p = [C(L(:,1),:) C(L(:,2),:)]'
  %
  warning('deform_mls is obsolete. Use deform.m, the class, instead.');
  % parse input
  V = varargin{1};
  F = varargin{2};
  C = varargin{3};
  mlsd = varargin{4};
  grid = false;

  ii = 5;
  while(ii <= size(varargin,2))
    if(strcmp(varargin{ii},'Lines'))
      ii = ii + 1;
      assert(ii<=size(varargin,2));
      L = varargin{ii};
      assert(strcmp(mlsd.constr,'lines'));
    elseif(strcmp(varargin{ii},'Image'))
      grid = true;
      ii = ii + 1;
      assert(ii<=size(varargin,2));
      im = varargin{ii};
      ii = ii + 1;
      assert(ii<=size(varargin,2));
      X = varargin{ii};
      ii = ii + 1;
      assert(ii<=size(varargin,2));
      Y = varargin{ii};
      grid_h = size(X,1);
      grid_w = size(X,2);
      %[X,Y] = meshgrid(1:grid_w,1:grid_h);
      %% list of points in grid
      V = [X(:) Y(:)];
    end
    ii = ii+1;
  end

  if(strcmp(mlsd.constr,'lines'))
    assert(1 == exist('L','var'));
    assert(size(L,1) == size(mlsd.p,2));
  end

  fprintf( ...
    ['\nCLICK a control point to visualize its corresponding weights ' ...
    'on the mesh.\n' ...
    'DRAG a control point to deform the mesh.\n\n']);


  %% make a new figure window
  %figure
  % clear current figure
  clf
  if(grid)
    gsh = warp( ...
      reshape(V(:,1),grid_h,grid_w), ...
      reshape(V(:,2),grid_h,grid_w), ...
      zeros(grid_h,grid_w), ...
      im);
  else
    % plot the original mesh
    tsh = trisurf(F,V(:,1),V(:,2),zeros(size(V,1),1), ...
      'FaceColor','interp');
  end
  % 2D view
  view(2);
  axis equal
  axis manual
  hold on;
  % plot the control points (use 3D plot and fake a depth offset by pushing
  % control points up in z-direction)
  C_plot = scatter3( ...
    C(:,1),C(:,2),0.1+0*C(:,1), ... 
    's','MarkerFaceColor','b', 'MarkerEdgeColor','b',...
    'LineWidth',2,'SizeData',100, ...
    'ButtonDownFcn',@oncontrolsdown);
  hold off;
  set(gca, 'visible', 'off');
  set(gcf,'Color','white');

  % keep track of window xmin, xmax, ymin, ymax
  win_min = min(V(:,1:2));
  win_max = max(V(:,1:2));
  % keep track of down position
  down_pos = [];
  % keep track of last two drag positions
  drag_pos = [];
  last_drag_pos = [];
  % keep track of mesh vertices at mouse down
  down_V = [];
  % keep track of index of selected control point
  ci = [];
  % keep track of control positions at mouse down
  global new_C;
  new_C = [];
  %% keep track of transformations stored at each control point, for 2D this is
  %% a 2x3 matrix
  %TR = repmat(eye(2,3),[1,1,size(C,1)]);
  % keep track of translations stored at each control point, for 2D this is a m
  % by 2 list of (x,y) vectors
  T = zeros(size(C,1),2);
  % keep track of rotations stored at each control point, for 2D this is a m
  % by 1 list of angles
  R = zeros(size(C,1),1);

  % type of click ('left','right')
  down_type  = '';

  function oncontrolsdown(src,ev)
    % get current mouse position
    down_pos=get(gca,'currentpoint');
    down_pos=[down_pos(1,1,1) down_pos(1,2,1)];
    last_drag_pos=down_pos;
    drag_pos=down_pos;
    % keep track of control point positions at mouse down
    new_C = [get(C_plot,'XData')' get(C_plot,'YData')'];
    % get index of closest control point
    [minD,ci] =  ...
      min(sum((new_C(:,1:2) - ...
      repmat(down_pos,size(new_C,1),1)).^2,2));
    % set color of mesh plot to weights of selected
    if(~strcmp(mlsd.constr,'lines'))
      if(~grid)
        set(tsh,'CData',mlsd.w(ci,:)');
      end
    end
    %% Set color to trinary visualization of 1-set, 0-set, and "inactive" set
    %WW = 0.5+zeros(size(W(:,ci)));
    %WW(abs(W(:,ci) - 1.0)<1e-1) = 1.0;
    %WW(abs(W(:,ci) - 0.0)<1e-1) = 0.0;
    %set(tsh,'CData',WW);
    % tell window that drag and up events should be handled by controls
    set(gcf,'windowbuttonmotionfcn',@oncontrolsdrag)
    set(gcf,'windowbuttonupfcn',@oncontrolsup)
    set(gcf,'KeyPressFcn',@onkeypress)
    if(strcmp('normal',get(gcf,'SelectionType')))
      % left-click
      down_type = 'left';
    else
      % other (right) click
      down_type = 'right';
    end
  end

  function oncontrolsdrag(src,ev)
    % keep last drag position
    last_drag_pos=drag_pos;
    % get current mouse position
    drag_pos=get(gca,'currentpoint');
    drag_pos=[drag_pos(1,1,1) drag_pos(1,2,1)];
    if(strcmp('left',down_type))
      % move selected control point by drag offset
      new_C(ci,:) = new_C(ci,:) + drag_pos-last_drag_pos;
      % update translation part of transformation stored at selected control
      T(ci,:) = T(ci,:) + (drag_pos-last_drag_pos);
    else
      R(ci) = R(ci) + 2*pi*(drag_pos(1)-last_drag_pos(1))/100;
    end

    % update mesh positions

    %% USING LINEAR BLEND SKINNING
    %TR = repmat(eye(2,3),[1,1,size(C,1)]);
    %TR(1,1,:) = cos(R);
    %TR(1,2,:) = -sin(R);
    %TR(2,1,:) = sin(R);
    %TR(2,2,:) = cos(R);
    %% stack TR as one tall 2*m by 2 matrix
    %RR = reshape(permute(TR(1:2,1:2,:),[2,1,3]),[2,2*size(C,1)])';
    %TR(1:2,3,:) = (T+C-stacktimes(RR,C))';
    %[new_V] = lbs(V,TR,W);

    %% USING DUAL QUATERNION SKINNING
    %% convert angles around z-axis to quaternions
    %Q = axisangle2quat(repmat([0,0,1],size(C,1),1),R);
    %RR = zeros(size(C,1)*2,2);
    %RR(1:2:end,1) =  cos(R);
    %RR(2:2:end,1) =  sin(R);
    %RR(1:2:end,2) = -sin(R);
    %RR(2:2:end,2) =  cos(R);
    %TT = [T + C-stacktimes(RR,C) zeros(size(T,1),1)];
    %% convert quaternions and translations into dualquaternions
    %DQ = quattrans2udq(Q,TT);
    %new_V = dualquatlbs(V,DQ,W);
    %new_V = new_V(:,1:2);

    update_positions(new_C);
  end

  function update_positions(new_C)
    % update display positions
    set(C_plot,'XData',new_C(:,1));
    set(C_plot,'YData',new_C(:,2));
    % MLS
    if(strcmp(mlsd.constr,'lines'))
      new_V = MLSD2DTransform(mlsd,[new_C(L(:,1),:) new_C(L(:,2),:)]')';
    else
      new_V = MLSD2DTransform(mlsd,new_C')';
    end


    if(grid)
      set(gsh,'XData',reshape(new_V(:,1),grid_h,grid_w));
      set(gsh,'YData',reshape(new_V(:,2),grid_h,grid_w));
    else
      set(tsh,'Vertices',new_V);
    end
  end

  function oncontrolsup(src,ev)
    % Tell window to handle drag and up events itself
    set(gcf,'windowbuttonmotionfcn','');
    set(gcf,'windowbuttonupfcn','');
    if(grid)
      X = get(gsh,'XData');
      Y = get(gsh,'YData');
      cur_V = [X(:) Y(:)];
    else
      cur_V = get(tsh,'Vertices');
    end
    % scale window to fit
    win_min = min([win_min; cur_V]);
    win_max = max([win_max; cur_V]);
    axis(reshape([win_min;win_max],1,2*size(cur_V,2)))
  end

  function onkeypress(src,ev)
    if(strcmp(ev.Character,'r'))
      update_positions(C);
    end
  end

end
