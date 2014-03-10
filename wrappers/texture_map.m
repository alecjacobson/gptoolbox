function [im_out,alpha_out] = texture_map(V,F,TC,im_in,alpha_in)
  % TEXTURE_MAP  texture map a give mesh (V,F) with a texture (TC,im_in) and save
  % it back to an image (im_out)
  %
  % [im_out,alpha_out] = texture_map(V,F,TC,im_in,alpha_in)
  %
  % Inputs:
  %   V  #V by 2 list of mesh vertex positions
  %   F  #F by 3 list of face indices
  %   TC  #V by 2 list of mesh texture coordinates, these should be between 0
  %     and 1, if they're not they will be normalized
  %   im_in  input texture image
  %   alpha_in  input texture image alpha channel
  % Outputs:
  %   im_out  deformed texture image
  %
  % Example:
  %  % original mesh (V,F) deformed mesh (NV,F)
  %  TC = texture_coords(V); 
  %  [oim, oalpha] = texture_map(V,F,TC,im,alpha);
  %

  if ~exist('alpha_in')
    alpha_in = ones(size(im_in,1),size(im_in,2));
  end
  if size(V,2) == 2
    V(:,end+1:3) = 0;
  end

  path_to_TexMapPreview = ...
    '/usr/local/igl/igl_apps/TexMapPreview/TexMapPreview';
  path_to_convert = '/opt/local/bin/convert';
  %path_to_TexMapPreviewConvert = ...
  %  '/usr/local/igl/igl_apps/TexMapPreview/TexMapPreviewConvert.sh';

  prefix = tempprefix();
  png_file_name = [prefix '.png'];
  tga_file_name = [prefix '.tga'];
  obj_file_name = [prefix '.obj'];
  output_tga_file_name = [prefix '-out.tga'];
  output_png_file_name = [prefix '-out.png'];

  % write image to png file
  imwrite(im_in,png_file_name,'Alpha',alpha_in);

  % convert the png image to a .tga image
  [status,result] = ...
    issue( ...
      'export DYLD_LIBRARY_PATH="";', ...
      path_to_convert, ...
      png_file_name, ...
      tga_file_name);
  % error on error
  if(status ~= 0)
      error(result);
  end
  % write mesh to obj file
  writeOBJ(obj_file_name,V,F,TC,F);

  % convert the png image to a .tga image
  [status,result] = ...
    issue( ...
      path_to_TexMapPreview, ...
      obj_file_name, ...
      tga_file_name, ...
      output_tga_file_name);
  % error on error
  if(status ~= 0)
      error(result);
  end

  % convert output tga to png

  % convert the png image to a .tga image
  [status, result] = ...
    issue( ...
      'export DYLD_LIBRARY_PATH="";', ...
      path_to_convert, ...
      output_tga_file_name, ...
      output_png_file_name);
  % error on error
  if(status ~= 0)
      error(result);
  end

  % clean up temporary files
  delete(tga_file_name);
  delete(png_file_name);
  delete(obj_file_name);
  delete(output_tga_file_name);

  if(nargout == 1)
    im_out = imread(output_png_file_name);
  elseif(nargout == 2)
    [im_out,m,alpha_out] = imread(output_png_file_name);
  else
    % clean up png file before error
    delete(output_png_file_name);
    error('Too many outputs');
  end
  % clean up png file
  delete(output_png_file_name);


  function [status, result] = issue(varargin)
  % ISSUE issue a system call with each arguement separated by a space
    command = [];
    for ii = 1:nargin
      command = ['' command ' ' varargin{ii}];
    end
    % show the command
    fprintf('%s\n',command);
    % issue the command
    [status, result] = system( command );
  end

end
