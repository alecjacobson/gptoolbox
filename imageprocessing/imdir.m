function files = imdir(dir_name)
  % IMDIR returns image files in directory
  %
  % files = imdir(dir_name)
  % 
  % Inputs:
  %   dir_name  path to directory
  % Outputs:
  %   files  struct array of image files in dir_name
  %
  % See also: dir
  %

  % get list of matlab readable file formats
  formats = imformats();
  % build list of formats' extensions: lower and upper case
  exts = [formats.ext upper([formats.ext])];

  files = [];
  for k = 1:size(exts,2);
    files = [files; dir([dir_name '/*.' exts{k}])];
  end
end
