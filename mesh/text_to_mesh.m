function [V,F,V0,F0] = text_to_mesh(str,varargin)
  % TEXT_TO_MESH Create an 3D mesh the given text
  %
  % [V,F] = text_to_mesh(str)
  % [V,F] = text_to_mesh(str,'ParameterName',ParameterValue, ...)
  % 
  % Inputs:
  %   str  string of text
  %   Optional:
  %     'Dilation'  followed by amount to dilate text image
  %     'FontName  followed by name of font. List available fonts using 
  %       `/usr/local/bin/convert -list font | sed -n "s/^ *Font: //p"` 
  %       {[]} (see `text_to_image`)
  %     'PointSize'  followed by height of font (in pixels) {200}
  % Outputs:
  %   V  #V by 3 list of mesh vertices
  %   F  #F by 3 list of face indices into V
  %
  % Example:
  %   im = text_to_image('Zapfino','FontName','Zapfino');imshow(im)
  % 

  % Defaults
  fontname = [];
  pointsize = 200;
  filename = 'text_to_image-tmp.png';
  r = 1;
  triangle_flags = '';
  siters=10;
  tol = 0.2;
  levels = 1;
  do_extrude = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Dilation','Extrude','FontName','PointSize','SmoothingIters','Tol','TriangleFlags'}, ...
    {'r','do_extrude','fontname','pointsize','siters','tol','triangle_flags'});
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

  im = text_to_image(str,'FontName',fontname,'PointSize',pointsize);
  [V0,F0] = bwmesh( ...
    imdilate(padarray(~im,[r,r]),strel('disk',r)), ...
    'Tol',tol, ...
    'SmoothingIters',siters, ...
    'TriangleFlags',triangle_flags);
  if do_extrude
    [V,F] = extrude(V0,F0,'Levels',levels);
  else
    V = V0;
    F = F0;
  end

end
