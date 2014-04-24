function G = depth_bars(im,D,d,varargin)
  % DEPTH_BARS Add vertical bars to an image at a specific depth.
  %
  % Inputs:
  %  im  h by w by c input image
  %  D  h by w 1./depth
  %  d  depth of bars
  %  Optional:
  %    'Color' followed by color of bars {0.8}
  % Outputs:
  %  G  h by w by c image with bars 
  %  

  % default values
  nc = size(im,3);
  color = repmat(0.8,[1 nc]);
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'Color'}, {'color'});
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
  color

  w = ceil(size(im,2)/40);
  B = im;
  for c = 1:nc
    B(:,[floor(end/3-w/2+(1:w)) floor(2*end/3-w/2+(1:w))],c) = color(c);
  end
  G = bsxfun(@times,B,D<d) + bsxfun(@times,im,D>=d);
end
