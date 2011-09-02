function [U,D] = unique_images(varargin)
  % UNIQUE_IMAGES 
  %
  % [U,D] = unique_images(directory_name)
  %
  % Inputs:
  %   directory_name
  % Outputs:
  %   U list of unique images/image paths
  %   D list of duplicate iamges, so that [U,D] are all images
  %

  dir_name = varargin{1};
  assert(isdir(dir_name));
  files = imdir(dir_name);
  U = {};
  ucount = 0;
  D = {};
  dcount = 0;
  valid = true(size(files));

  for ii = 1:numel(files)
    fii = [dir_name '/' files(ii).name];
    imii = imread(fii);
    for jj = (ii+1):numel(files)
      imjj = imread([dir_name '/' files(jj).name]);
      if all(size(imii) == size(imjj)) && all(imii(:) == imjj(:))
        valid(jj) = false;
      end
    end
  end
  U = files(valid);
  D = files(~valid);

end
