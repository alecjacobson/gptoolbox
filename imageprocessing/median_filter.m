function M = median_filter(im)
  % MEDIAN_FILTER simpler wrapper for calling medfilt2 on each channel
  %
  % M = median_filter(im)
  %
  % Input:
  %   im  w by h by c image
  % Output:
  %   M  median filtered image in each channel
  %

  % This could probably be a one-liner using num2cell and cellfun
  M = zeros(size(im));
  for c = 1 : size(im,3)
    M(:,:,c) = medfilt2(im(:,:,c));
  end
end
