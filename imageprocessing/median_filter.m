function M = median_filter(im)
  % MEDIAN_FILTER simpler wrapper for calling medfilt2 on each channel
  % M = median_filter(im)
  M = zeros(size(im));
  for c = 1 : size(im,3)
    M(:,:,c) = medfilt2(im(:,:,c));
  end
end
