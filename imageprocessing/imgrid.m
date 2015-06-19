function varargout = imgrid(im,varargin)
  % IMGRID Display an image with a overlay grid.
  %
  % t = imgrid(im)
  %
  % Inputs:
  %   im  h by w by channels list of colors/values
  % Outputs:
  %   t  handle to surface with colored faces
  %

  corner = [1 1];
  show_values = false;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'ShowValues'},{'show_values'});
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

  [X,Y] = meshgrid(-0.5+(1:(size(im,2)+1)),-0.5+(1:(size(im,1)+1)));
  s = surf(X,Y,zeros(size(X)));
  [Q,V] = surf2patch(s);
  if size(im,3) == 1
    t = trisurf(Q,V(:,1),V(:,2),V(:,3),'CData',im(:),'EdgeColor',[0.2 0.2 0.2]);
    if show_values
      [BX,BY] = meshgrid(1:(size(im,2)),1:(size(im,1)));
      nz = find(im);
      text(BX(nz),BY(nz),num2str(im(nz)), ...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',16);
    end
  else
    imC = reshape(im,[],3);
    t = trisurf(Q,V(:,1),V(:,2),V(:,3),'FaceVertexCData',imC);
  end
  axis equal;
  set(gca,'YDir','reverse');
  set(gca,'XAxisLocation','top', ...
    'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4]);
  view(2);
  if nargout>0
    varargout{1} = t;
  end
end
