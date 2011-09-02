function imshow_rgb(R,G,B)
  % Wrapper for imshow taking R,G,B values separately


  % compose the r,g,b signals into 3d matrix
  rgb = ziprgb(R,G,B);
  % show the image
  imshow(rgb);

  clear rgb
end
