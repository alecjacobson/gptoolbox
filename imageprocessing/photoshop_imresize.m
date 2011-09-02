function rim = photoshop_imresize(im,z)
  % PHOTOSHOP_IMRESIZE resize an image like photoshop with a given zoom. If
  % zoomed in enough then faint lines appear separating pixels
  %
  % rim = photoshop_imresize(im,z)
  %
  % Inputs:
  %   im  image to show
  %   z  zoom amount, 
  % Outputs:
  %   rim  resized with pixel separators if necesary
  %

  % zoom scale should be either smaller than one or an integer
  if z > 1
    z = round(z);
  end

  if z <= 5
    rim = imresize(im,z);
  else
    % number of channels
    nc = size(im,3);
    % resize image according to zoom scale with nearest neighbor method
    rim = imresize(im,z,'nearest');
    % add faint pixel separators, backwards engineered from photoshop
    % build mask
    pixel_sep_mask = uint8(zeros(z,z));
    pixel_sep_mask(1,1) = 60;
    pixel_sep_mask(1,end) = 65;
    pixel_sep_mask(end,end) = 60;
    pixel_sep_mask(end,1) = 55;
    pixel_sep_mask(1,2:(end-1)) = 35;
    pixel_sep_mask(end,2:(end-1)) = 29;
    pixel_sep_mask(2:(end-1),1) = 29;
    pixel_sep_mask(2:(end-1),end) = 35;

    if isfloat(rim)
      % convert mask to double
      pixel_sep_mask = im2double(pixel_sep_mask);
      one = 1.0;
    elseif islogical(rim)
      % convert mask to double
      pixel_sep_mask = im2double(pixel_sep_mask);
      % convert image to double
      rim = 1.0*rim;
      one = 1.0;
    else
      one = 255;
    end

    pixel_sep_mask = repmat(pixel_sep_mask,size(im));

    if strcmp('uint8',class(rim));
      rim = uint8(uint16(rim).*uint16((one-pixel_sep_mask))/one) + pixel_sep_mask;
    else
      rim = rim.*one-pixel_sep_mask + pixel_sep_mask;
    end

  end

end
