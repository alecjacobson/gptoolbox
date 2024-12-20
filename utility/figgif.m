function f = figgif(filename,varargin)
  % FIGGIF Write the current figure to filename. If filename exists, append it
  % as an animation frame with zero delay
  %
  % Inputs:
  %   filename  path to gif file
  % Outputs:
  %   f  flag whether file already existed
  % 

  if isempty(filename)
    f = [];
    return;
  end

  rgb2ind_args = {};
  delay = 0;
  cdata = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'CData','Delay','rgb2ind'}, ...
    {'cdata','delay','rgb2ind_args'});
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

  if isempty(cdata)
    frame = getframe(gcf);
    cdata = frame.cdata;
  else
    if size(cdata,3) == 1
      cdata = cat(3,cdata,cdata,cdata);
    end
    cdata = im2uint8(cdata);
  end
  [SIf,cm] = rgb2ind(cdata,256,rgb2ind_args{:});
  f = exist(filename,'file');
  if ~f
    imwrite(SIf,cm,filename,'Loop',Inf,'Delay',delay);
  else
    imwrite(SIf,cm, filename,'WriteMode','append','Delay',delay);
  end
end
