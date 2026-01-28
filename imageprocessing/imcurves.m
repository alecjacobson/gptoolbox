function [Z,XV] = imcurves(im,XV,varargin)
  % IMCURVES  Filter an image using a spline curve similar to Photoshop's curves
  % tool
  %
  % [Z,XV] = imcurves(im,XV)
  % [Z,XV] = imcurves(im,XV,'ParameterName',ParameterValue,...)
  % 
  % Inputs:
  %   im  h by w by c image
  %   XV  n by 2 list of spline parameter locations and value pairs
  %   Optional:
  %     'Interactive'  followed by true or false, if true then the user can
  %       interactively add and move control points
  % Outputs:
  %   Z  h by w by c image after filtering
  %   XV  n by 2 list of spline parameter locations and value pairs
  interactive = true;
  if nargin<2 || isempty(XV)
    XV = [[0;1] [0;1]];
    interactive = true;
  else
    interactive = false;
  end


  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Interactive'}, ...
    {'interactive'});
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

  if interactive
    pos = get(gcf,'Position');
    if pos(3)*2<size(im,2)
      warning('interactive display on smaller image');
      im_fit = imresize(im,pos(3)/size(im,2));
    else
      im_fit = im;
    end
    Z_fit = filter(im_fit);
    subplot(1,2,1);
    ish = imshow(Z_fit);
    subplot(1,2,2);
    t = linspace(0,1,100)';
    c = filter(t);
    psh = plot(t,c,'LineWidth',3);
    hold on;
    ssh = sct(XV,'w','filled','MarkerEdgeColor','k');
    hold off;
    axis equal;
    axis(0.5+1.25*[-0.5 0.5 -0.5 0.5]);



    % Add callbacks so that if user clicks on ssh then the corresponding row in XV
    % is moved and the plot updated.
    set(gcf,'WindowButtonMotionFcn','');
    set(ssh,'ButtonDownFcn',@startDragFcn);
    set(psh,'ButtonDownFcn',@startNewDragFcn);
    set(gcf,'WindowButtonUpFcn',@stopDragFcn);
    selected_index = [];

    % Add callback so if user presses delete the selected index is deleted
    set(gcf,'KeyPressFcn',@keyPressFcn);

    % wait for state.finished to be true
    set(ssh,'UserData','waiting')
    waitfor(ssh,'Userdata','finished');
  else
    Z = filter(im);
  end

  function finish()
    set(gcf,'KeyPressFcn','');
    set(gcf,'WindowButtonMotionFcn','');
    set(ssh,'ButtonDownFcn','');
    set(gcf,'WindowButtonUpFcn','');
    Z = filter(im);
  end

  function keyPressFcn(src,event)
    if strcmp(event.Key,'return')
      finish();
      set(ssh,'UserData','finished');
    end
    if strcmp(event.Key,'backspace')
      if ~isempty(selected_index) && size(XV,1)>2
        XV(selected_index,:) = [];
        selected_index = [];
        update_plot();
      end
    end
  end

  function Y = filter(Z)
    Y = interp1(XV(:,1),XV(:,2),Z,'spline');
    [~,xi] = max(XV(:,1));
    [~,ni] = min(XV(:,1));
    Y(Z>XV(xi,1)) = XV(xi,2);
    Y(Z<XV(ni,1)) = XV(ni,2);
    Y = max(min(Y,1),0);
  end

  function update_plot()
    c = filter(t);
    ish.CData = filter(im_fit);
    set(psh,'YData',c);
    set(ssh,'XData',XV(:,1),'YData',XV(:,2));
  end
  function startNewDragFcn(varargin)
    set(gcf,'WindowButtonMotionFcn',@draggingFcn);
    cp = get(gca,'CurrentPoint');
    cp = cp(1,1:2);
    index_below = find(XV(:,1)<cp(1),1,'last');
    XV = [XV(1:index_below,:); cp; XV(index_below+1:end,:)];
    selected_index = index_below+1;
    update_plot();
  end
  function startDragFcn(varargin)
    set(gcf,'WindowButtonMotionFcn',@draggingFcn);
    cp = get(gca,'CurrentPoint');
    cp = cp(1,1:2);
    [~,selected_index] = min(sum((XV-cp).^2,2));
  end
  function draggingFcn(varargin)
    if isempty(selected_index)
      return;
    end
    cp = get(gca,'CurrentPoint');
    cp = max(min(cp(1,1:2),1),0);
    XV(selected_index,:) = cp;
    update_plot();
  end
  function stopDragFcn(varargin)
    set(gcf,'WindowButtonMotionFcn','');
  end
end
