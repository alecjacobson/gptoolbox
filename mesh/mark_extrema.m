function [s] = mark_extrema(varargin)
  % MARK_EXTREMA mark with colored points the extrema of the colordata on a
  % given surface plot
  %
  % [s] = mark_extrema(f,'ParameterName',ParameterValue);
  % [s] = mark_extrema(f,S,'ParameterName',ParameterValue);
  %
  % Inputs:
  %   t  handle to trisurf plot
  %   Optional:
  %     S  optional scalar data
  %     'Epsilon' followed by epsilon value used for near extrema
  %     'Zero' followed by zero value used for extrema
  %     'MaxNum' maximum number of extrema to show
  % Outputs:
  %   s  list of handles to scatter3 plots for extrema
  %
  t = varargin{1};
  ii = 2;
  if nargin > 1 && ~ischar(varargin{2})
    S = varargin{ii};
    ii = ii + 1;
  end

  % defaults
  epsilon = 1e-15;
  max_num = Inf;
  my_zero = 0;

  while(ii <= nargin)
    switch varargin{ii}
    case 'Zero'
      ii = ii + 1;
      assert(ii<=nargin);
      my_zero = varargin{ii};
    case 'Epsilon'
      ii = ii + 1;
      assert(ii<=nargin);
      epsilon = varargin{ii};
    case 'MaxNum'
      ii = ii + 1;
      assert(ii<=nargin);
      max_num = varargin{ii};
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii+1;
  end


  if strcmp(get(t,'type'),'line')
    V = [get(t,'XData');get(t,'YData');0*get(t,'XData')]';
    F = [1:(size(V,1)-1); 2:size(V,1)]';
    SD = get(t,'YData')';
  else
    V = get(t,'Vertices');
    F = get(t,'Faces');
    SD = get(t,'FaceVertexCData');
  end

  if ~exist('S','var')
    S = SD;
  end

  [xv,xi,xd] = local_max(F,S);
  [nv,ni,nd] = local_min(F,S);
  xi = xi(abs(xd)>my_zero);
  ni = ni(abs(nd)>my_zero);

  if max_num ~= Inf
    xi = unique([xi(xv==max(xv));xi(1:round(max(1,end/max_num)):end)]);
    ni = unique([ni(nv==min(nv));ni(1:round(max(1,end/max_num)):end)]);
  end

  [~,wxi,wxd] = local_max(F,S,epsilon);
  [~,wni,wnd] = local_min(F,S,epsilon);
  wxi = wxi(abs(wxd)>my_zero);
  wni = wni(abs(wnd)>my_zero);

  wxi = setdiff(wxi,xi);
  wni = setdiff(wni,ni);

  flati = intersect([xi;wxi],[ni;wni]);
  xi = setdiff(xi,flati);
  ni = setdiff(ni,flati);
  wxi = setdiff(wxi,flati);
  wni = setdiff(wni,flati);

  size_data = 200;

  show_flat = false;

  a = get(t,'Parent');
  hold on;
  s = [];
  s = [s(:);point(a,V(xi,:),[16 2 1]/16,size_data)];
  s = [s(:);point(a,V(wxi,:),[16 10 10]/16,size_data)];
  s = [s(:);point(a,V(ni,:),[1 2 16]/16,size_data)];
  s = [s(:);point(a,V(wni,:),[10 10 16]/16,size_data)];
  if show_flat
    s = [s(:);point(a,V(flati,:),[16 16 10]/16,size_data)];
  end
  hold off;
end
