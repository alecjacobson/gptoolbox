function f = figpng(filename,varargin)
  % FIGPNG Write the current figure to filename. If filename exists, replace it
  %
  % Inputs:
  %   filename  path to png file
  % Outputs:
  %   f  flag whether file already existed
  % 

  alpha = false; 
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Alpha'}, ...
    {'alpha'});
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




  if alpha
    old_bg = get(gcf,'Color');

    set(gcf,'Color','w');
    IW = im2double(getfield(getframe(gcf),'cdata'));
    set(gcf,'Color','k');
    IK = im2double(getfield(getframe(gcf),'cdata'));
    A = mean((IK - IW) + 1,3);
    I = IK./A;

    % reset background color
    set(gcf,'Color',old_bg);

    imwrite(I,filename,'Alpha',A);
  else
    frame = getframe(gcf);
    imwrite(frame.cdata,filename);
  end
end

