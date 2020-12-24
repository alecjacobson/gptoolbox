function rgb = oklab2rgb(ok)
  % OKLAB2RGB Convert colors stored in OKLab color space to RGB color space
  %
  % rgb = oklab2rgb(ok)
  %
  % Inputs:
  %   ok  #h by #w by 3 image of colors or #ok by 3 list of colors
  % Outputs:
  %   rgb  converted ok colors (same size as ok)
  %
  % See also: rgb2lin, lin2rgb, oklab2lin , rgb2oklab, lin2oklab
  lin = oklab2lin(ok);
  rgb = lin2rgb(lin);
end

