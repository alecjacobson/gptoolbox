function im = text_to_image(str,varargin)
  % TEXT_TO_IMAGE  Create an grayscale image the given text
  %
  % im = text_to_image(str)
  % im = text_to_image(str,'ParameterName',ParameterValue, ...)
  % 
  % Inputs:
  %   str  string of text
  %   Optional:
  %     'FontName  followed by name of font. List available fonts using 
  %       `/usr/local/bin/convert -list font | sed -n "s/^ *Font: //p"` 
  %       {'CourierNewB'}
  %     'PointSize'  followed by height of font (in pixels) {200}
  % Outputs:
  %   im  h by w image of text
  %
  % Example:
  %   im = text_to_image('Zapfino','FontName','Zapfino');imshow(im)
  % 

  % Defaults
  fontname = [];
  pointsize = 200;
  filename = 'text_to_image-tmp.png';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'FontName','PointSize'}, ...
    {'fontname','pointsize'});
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
  if isempty(fontname)
    fontname = 'CourierNewB';
  end

  cmd = sprintf( ...
    ['%s -size %dx%d xc:white -font "%s" -pointsize %d ' ...
     '-fill black -annotate +%d+%d "%s" -trim -bordercolor "#FFF" ' ...
     '-border 0 %s'], ...
   path_to_convert,pointsize*2*numel(str),pointsize*4,fontname,pointsize, ...
   0,2*pointsize,str, ...
   filename);
  [status,result] = system(cmd);
  if status ~= 0
    warning(cmd);
    error(result);
  end
  im = imread(filename);
  delete(filename);
end
