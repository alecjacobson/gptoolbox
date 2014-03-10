function [U,D] = unique_images(dir_name)
  % UNIQUE_IMAGES Generate a list of unique and duplicate images in a
  % directory.
  %
  % [U,D] = unique_images(dir_name)
  %
  % Inputs:
  %   dir_name  path to directory
  % Outputs:
  %   U list of unique images/image paths
  %   D list of duplicate iamges, so that [U,D] are all images
  %

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
