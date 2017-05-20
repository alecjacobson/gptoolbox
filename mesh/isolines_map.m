function M = isolines_map(varargin)
  % ISOLINES_MAP Create a colormap that can be used to create "isolines". That
  % is colored intervals where everyother interval is a thin black interval.
  % 
  % M = isolines_map(n)
  % M = isolines_map(N)
  % M = isolines_map(...,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   n  number of *isolines* (thus there will be n+1 colored intervals), using
  %   default colormap
  %   or
  %   N  #N by 3 input colormap, (thus there will be #N-1 isolines)
  %   Optional:
  %     'LineColor'  followed by {1|#N-1} by 3 color of isolines {[0 0 0]}
  %     'Thickness'  ratio of isointerval thick to isoline, where thick is
  %       measured on the color axis {10}
  %
  % Note: Remember to use 'FaceLighting','phong' in your trisurf plot to get
  %   crisp isointervals (and isolines)
  %

  if numel(varargin{1}) == 1
    n = varargin{1};
    N = parula(n+1);
  else
    N = varargin{1};
    assert(size(N,2) == 3,'Number columns in colormap (%d) not 3',size(N,2));
    n = size(N,1) - 1;
  end

  % Default optional parameters
  thick = 10;
  line_color = zeros([n 3]);

  ii = 2;
  while ii <= nargin
    switch varargin{ii}
    case 'Thickness'
      ii = ii + 1;
      assert(ii<=nargin);
      thick = varargin{ii};
      assert(thick == round(thick),'Thickness (%d) ot an integer',thick);
    case 'LineColor'
      ii = ii + 1;
      assert(ii<=nargin);
      line_color = varargin{ii};
      assert(size(line_color,2) == 3, ...
        'Number of columns in LineColor %d not 3',size(line_color,2));
      if size(line_color,1) == 1
        line_color = repmat(line_color,[n 1]);
      else
        assert(size(line_color,1) == n, ...
          'Number of rows in LineColor (%d) not n (%d)',size(line_color,1),n);
      end
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii+1;
  end

  % append phony color to line_color to make it match N
  line_color = [line_color;[1 0 0]];
  assert(size(N,1) == size(line_color,1))

  % repeat N thick times
  N = repmat(N,[1 1 thick]);
  M = reshape(permute(cat(3,N,line_color),[3 1 2]),(n+1)*(thick + 1),3);
  % drop phony color 
  M = M(1:(end-1),:);

end
