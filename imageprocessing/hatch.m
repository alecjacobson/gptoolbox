function H = hatch(im,levels,thickness)
  % HATCH computes dithering hatched image with a number of levels of specified
  % thickness
  %
  % H = hatch(im,levels)
  %
  % Inputs:
  %  im  original image
  %  levels  number of hatching levels {3}
  %  thickness  thickness of individual hatch lines
  % Outputs:
  %  H  hatched image
  %
  % See also: dith
  %

  assert(levels>1);
  
  % make sure image is grayscale
  switch size(im,3)
  case 1
  % do nothing, already "grayscale"
  case 3
    warning('Converting input image to grayscale using rgb2gray');
    im = rgb2gray(im);
  otherwise
    error('Input image should be grayscale');
  end

  % make sure image is double
  if ~isfloat(im)
    warning([ ...
      'Converting image to double, using im2double, ' ...
      'output will also be double']);
    im = im2double(im);
  end


  H = zeros(size(im));
  for i = 0:(levels-1)
    % compute ratio of black to white for this level
    ratio = i/(levels-1);
    lo = i/levels;
    % compile pattern with this ratio of black to white
    pattern = zeros(1,levels-1)';
    for j = 1:(levels-1)
      pattern(j) = sum(pattern(1:j))/j < ratio;
    end
    pattern = repmat(pattern,1,thickness)';
    pattern = pattern(:);
    scale = ceil(size(im)./size(pattern));
    pattern_mask = repmat(pattern,scale(1),scale(2));
    pattern_mask = pattern_mask(1:size(im,1),1:size(im,2));
    %if(i == round((levels-1)/2))
    %  imshow(pattern_mask);
    %end
    H( im > lo ) = pattern_mask(im > lo);
  end
  % get rid of small lines 
  % b = 1-bwareaopen(1-bwareaopen(h,10),10);
  % get rid of speckles in lines
  % m = medfilt2(b,[1,5]);
end
