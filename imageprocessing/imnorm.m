function N = imnorm(im)
  % IMNORM
  % 
  % N = imnorm(im)
  %
  % Normalize image to be between the range 0 and 1. Current just works with
  % double images.
  %
  % Inputs:
  %   im  original input image
  % Output:
  %   N  normalized image
  %
  %
  warning( [...
     'This is a legacy wrapper... ' ...
     'You should call matrixnormalize directly.\n']);
  N = matrixnormalize(im);
end
