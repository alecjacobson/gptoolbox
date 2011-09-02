% RES = blur(IM, LEVELS, FILT)
%
% Blur an image, by filtering and downsampling LEVELS times
% (default=1), followed by upsampling and filtering LEVELS times.  The
% blurring is done with filter kernel specified by FILT (default =
% 'binom5'), which can be a string (to be passed to namedFilter), a
% vector (applied separably as a 1D convolution kernel in X and Y), or
% a matrix (applied as a 2D convolution kernel).  The downsampling is
% always by 2 in each direction.

% Eero Simoncelli, 3/04.

function res = blur(im, nlevs, filt)

%------------------------------------------------------------
%% OPTIONAL ARGS:

if (exist('nlevs') ~= 1) 
  nlevs = 1;
end

if (exist('filt') ~= 1) 
  filt = 'binom5';
end

%------------------------------------------------------------

if isstr(filt)
  filt = namedFilter(filt);
end  

filt = filt/sum(filt(:));

if nlevs > 0
  if (any(size(im)==1))
    if (~any(size(filt)==1))
      error('Cant  apply 2D filter to 1D signal');
    end
    if (size(im,2)==1)
      filt = filt(:);
    else
      filt = filt(:)';
    end
    
    in = corrDn(im,filt,'reflect1',(size(im)~=1)+1);
    out = blur(in, nlevs-1, filt);
    res = upConv(out, filt, 'reflect1', (size(im)~=1)+1, [1 1], size(im));

  elseif (any(size(filt)==1))
    filt = filt(:);

    in = corrDn(im,filt,'reflect1',[2 1]);
    in = corrDn(in,filt','reflect1',[1 2]);
    out = blur(in, nlevs-1, filt);
    res = upConv(out, filt', 'reflect1', [1 2], [1 1], [size(out,1),size(im,2)]);
    res = upConv(res, filt, 'reflect1', [2 1], [1 1], size(im));

  else

    in = corrDn(im,filt,'reflect1',[2 2]);
    out = blur(in, nlevs-1, filt);
    res = upConv(out, filt, 'reflect1', [2 2], [1 1], size(im));
  end
else
  res = im;
end

