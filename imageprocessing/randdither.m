function D = randdither(im)
  % RANDDITHER Dither an image using random noise
  %
  % D = randdither(im)
  %
  % Input:
  %   im  grayscale image.
  % Output:
  %   D  dithered image
  %
  % See also: dith


  h = size(im,1);
  w = size(im,2);

  K = [1,2,1;2,0,2;1,2,1];
  offset = [1,1];
  K = fspecial('gaussian',5,5);
  offset = [2,2];

  K(1+offset(1),1+offset(2)) = 0;
  K = K./sum(K(:));

  D = zeros(size(im)+size(K));
  O = zeros(size(im)+size(K));

  % original image
  O(offset(1)+(1:h),offset(2)+(1:w)) = im;

  % random path through pixels of original image
  [I,J] = ind2sub([h,w],1:(h*w)); 
  P = [I;J]';
  P = P(randperm(size(P,1)),:);
  P(:,1) = P(:,1) + offset(1);
  P(:,2) = P(:,2) + offset(2);

  % keep track of which pixels have not been processed
  F = ones(size(im)+size(K));

  % precompute indices of diffusion region
  dY = (1:size(K,1))-1-offset(1);
  dX = (1:size(K,2))-1-offset(2);

  for pi = 1:size(P,1)
    y = P(pi,1);
    x = P(pi,2);
    D(y,x) = round(O(y,x));
    e = O(y,x) - D(y,x);
    F(y,x) = 0;
    Y = y + dY;
    X = x + dX;
    KK = F(Y,X).*K;
    if(sum(KK(:))>0)
      % commenting this normalization out, has sort of an Atkinson effect
      KK = KK./sum(KK(:));
      O(Y,X) = O(Y,X) + KK.*e;
    end
  end

  % resize D to be same as original image
  D = D(offset(1)+(1:h),offset(2)+(1:w));

end
