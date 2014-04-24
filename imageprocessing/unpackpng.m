function [G,depth] = unpackpng(varargin)
  % UNPACKPNG  unpack a rgba image read from a png file to a grayscale image of
  % doubles
  %
  % [G,depth] = unpackpng(filename)
  % [G,depth] = unpackpng(rgba)
  %
  % Inputs:
  %   filename  path to png file
  %   rgba  h by w by 4  RGBA image
  % Outputs:
  %   G  h by w image of doubles
  %   depth  number of bits used: 8 or 16
  %

  if ischar(varargin{1})
    [rgb,~,alpha] = imread(varargin{1});
    rgba = cat(3,rgb,alpha);
  else
    rgba = varargin{1};
    if ~isa(rgba,'uint16')
      warning('Converting input to uint16');
      rgba = im2uint16(rgba);
    end
  end

  if isa(rgba,'uint16')
    depth = 16;
  elseif isa(rgba,'uint8')
    depth = 8;
  else
    error('Unable to recognize bit depth');
  end
  %rgba = double(rgba);
  %rgba = im2double(rgba)*(1-2^-depth)*2^depth;
  rgba = im2double(rgba);
  squeeze(rgba(1,:,:))
  rgba = rgba*(1-2^-depth)*2^depth;
  % unpack
  G = zeros(size(rgba,1),size(rgba,2));
  for ii = 1:4
    G = G + rgba(:,:,ii)*2^(-ii*depth);
  end
  G = G/(1 - 2^-depth);

end
