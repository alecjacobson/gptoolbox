function N = matrixnormalize(M)
  % MATRIXNORMALIZE Normalize matrix values to be between the range 0 and 1.
  % Current just works with matrices of type double.
  % 
  % N = matrixnormalize(M)
  %
  % Inputs:
  %   M  original input matrix
  % Output:
  %   N  normalized matrix
  % 
  % Example:
  %   imshow(matrixnormalize(im))
  %   % Equivalent to
  %   imshow(im,[])
  %
  switch class(M)
  case {'double','single'}
    N = (M-min(M(:)))./(max(M(:))-min(M(:)));
  case 'logical'
    N = M;
  case 'uint8'
    N = (M-min(M(:)))*(255/double((max(M(:))-min(M(:)))));
  otherwise
    error('Class not supported');
  end
end
