function ok = rgb2oklab(rgb)
  % RGB2OKLAB Convert colors stored in RGB color space to OKLab color space
  %
  % ok = rgb2oklab(rgb)
  %
  % Inputs:
  %   rgb  #h by #w by 3 image of colors or #rgb by 3 list of colors
  % Outputs:
  %   ok  converted ok colors (same size as rgb)
  %
  % See also: rgb2lin, lin2rgb, oklab2rgb, lin2oklab, oklab2lin
  %
  % Example:
  %   C = permute([0.99215,0.99607,0.99607;0.14901,0.05490,0.94901],[1 3 2]);
  %   L = rgb2oklab(C);
  %   imshow([ ...
  %    interp1([0 1],C,linspace(0,1,334)).*ones(1,146) ...
  %    oklab2rgb(interp1([0 1],L,linspace(0,1,334)).*ones(1,146))]);
  lin = rgb2lin(rgb);
  ok = lin2oklab(lin);
end
