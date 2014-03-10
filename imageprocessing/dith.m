function D = dith(im,K,offset)
  % DITH Dithers an input image
  %
  % D = dither(im,K) Image im is dithered using algorithm specified by K:
  %   ['Floyd-Steinberg','Jarvis','Stucki','Sierra3','Sierra2','Sierra-lite',
  %    'Atkinson']
  %
  % D = dither(im,K,offset) Image im is dithered using provided error diffusion
  %   kernel K where current pixel is located in K by (1,1) + offset
  %
  % D = dither(im)  Image im is dithered using default algorithm
  %   {'Floyd-Steinberg'}
  %
  % Inputs:
  %   im  grayscale image to be dithered
  %   K  algorithm or kernel
  %   offset  pixel offset for kernel
  % Output:
  %   D  dithered image
  %

  % get width and height
  h = size(im,1);
  w = size(im,2);

  % Error diffusion method
  % error diffusion kernel and
  % offset from (1,1) of current pixel in kernel

  if(~exist('K','var'))
    K = 'Floyd-Steinberg';
  end
  
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

  if(strcmp(class(K),'char'))
    if(strcmp(K,'Floyd-Steinberg'))
      offset = [0,1];
      K = [0,0,7;3,5,1]/16;
    elseif(strcmp(K,'Jarvis'))
      offset = [0,2];
      K = [0 0 0 7 5;3 5 7 5 3;1 3 5 3 1]/48;
    elseif(strcmp(K,'Stucki'))
      offset = [0,2];
      K = [0 0 0 8 4; 2 4 8 4 2;1 2 4 2 1]/42;
    elseif(strcmp(K,'Sierra3'))
      offset = [0,2];
      K = [0 0 0 5 3;2 4 5 4 2;0 2 3 2 0]/32;
    elseif(strcmp(K,'Sierra2'))
      offset = [0,2];
      K = [0 0 0 4 3;1 2 3 2 1]/16;
    elseif(strcmp(K,'Sierra-lite'))
      offset = [0,0];
      K = [0,2;1,1]/4;
    elseif(strcmp(K,'Atkinson'))
      offset = [0,1];
      K = [0 0 1 1;1 1 1 0;0 1 0 0]/8;
    else
      error 'Unsupported algorithm name'
    end
  else
    assert(exist('offset','var')==1);
  end

  D = zeros(size(im)+size(K));
  O = zeros(size(im)+size(K));
  % original image
  O(offset(1)+(1:h),offset(2)+(1:w)) = im;

  % precompute indices of diffusion region
  dY = (1:size(K,1))-1-offset(1);
  dX = (1:size(K,2))-1-offset(2);

  % loop from top to bottom
  for y = (offset(1) + (1:h))
    Y = y+dY;
    % loop from left to right
    for x = (offset(2) + (1:w))
      D(y,x) = round(O(y,x));
      % get error at this pixel
      e = O(y,x) - D(y,x);
      % distribute error to nearby pixels using kernel
      X = x+dX;
      % diffuse error
      O(Y,X) = O(Y,X) + K.*e;
    end
  end

  % resize D to be same as original image
  D = D(offset(1)+(1:h),offset(2)+(1:w));

end
