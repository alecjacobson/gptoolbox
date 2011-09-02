function [mssim, ssim_map] = ssim_index(img1, img2, K, window, L)

%========================================================================
%Edited code by Adam Turcotte and Nicolas Robidoux
%Laurentian University
%Sudbury, ON, Canada
%Last Modified: 2011-01-22
%----------------------------------------------------------------------
%This code implements a refactored computation of SSIM that requires
%one fewer blur (4 instead of 5), the same number of pixel-by-pixel
%binary operations (10), and two fewer unary operations (6 instead of 8).
%As a result, it is about 20% faster.
%
%In addition, this version reduces memory usage with in-place functions.
%As a result, it supports larger input images.
%========================================================================

%========================================================================
%SSIM Index, Version 1.0
%Copyright(c) 2003 Zhou Wang
%All Rights Reserved.
%
%The author is with Howard Hughes Medical Institute, and Laboratory
%for Computational Vision at Center for Neural Science and Courant
%Institute of Mathematical Sciences, New York University.
%
%----------------------------------------------------------------------
%Permission to use, copy, or modify this software and its documentation
%for educational and research purposes only and without fee is hereby
%granted, provided that this copyright notice and the original authors'
%names appear on all copies and supporting documentation. This program
%shall not be used, rewritten, or adapted as the basis of a commercial
%software or hardware product without first obtaining permission of the
%authors. The authors make no representations about the suitability of
%this software for any purpose. It is provided "as is" without express
%or implied warranty.
%----------------------------------------------------------------------
%
%This is an implementation of the algorithm for calculating the
%Structural SIMilarity (SSIM) index between two images. Please refer
%to the following paper:
%
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error measurement to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 1, Jan. 2004.
%
%Kindly report any suggestions or corrections to zhouwang@ieee.org
%
%----------------------------------------------------------------------
%
%Input : (1) img1: the first image being compared
%        (2) img2: the second image being compared
%        (3) K: constants in the SSIM index formula (see the above
%            reference). defualt value: K = [0.01 0.03]
%        (4) window: local window for statistics (see the above
%            reference). default widnow is Gaussian given by
%            window = fspecial('gaussian', 11, 1.5);
%        (5) L: dynamic range of the images. default: L = 255
%
%Output: (1) mssim: the mean SSIM index value between 2 images.
%            If one of the images being compared is regarded as 
%            perfect quality, then mssim can be considered as the
%            quality measure of the other image.
%            If img1 = img2, then mssim = 1.
%        (2) ssim_map: the SSIM index map of the test image. The map
%            has a smaller size than the input images. The actual size:
%            size(img1) - size(window) + 1.
%
%Default Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   [mssim ssim_map] = ssim_index(img1, img2);
%
%Advanced Usage:
%   User defined parameters. For example
%
%   K = [0.05 0.05];
%   window = ones(8);
%   L = 100;
%   [mssim ssim_map] = ssim_index(img1, img2, K, window, L);
%
%See the results:
%
%   mssim                        %Gives the mssim value
%   imshow(max(0, ssim_map).^4)  %Shows the SSIM index map
%
%========================================================================


if (nargin < 2 || nargin > 5)
   ssim_index = -Inf;
   ssim_map = -Inf;
   return;
end

if (size(img1) ~= size(img2))
   ssim_index = -Inf;
   ssim_map = -Inf;
   return;
end

[M N] = size(img1);

if (nargin == 2)
   if ((M < 11) || (N < 11))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);	%
   K(1) = 0.01;                             % default settings
   K(2) = 0.03;                             %
   L = 255;                                 %
end

if (nargin == 3)
   if ((M < 11) || (N < 11))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 || K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 4)
   [H W] = size(window);
   if ((H*W) < 4 || (H > M) || (W > N))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 || K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 5)
   [H W] = size(window);
   if ((H*W) < 4 || (H > M) || (W > N))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   if (length(K) == 2)
      if (K(1) < 0 || K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));
img1 = double(img1);
img2 = double(img2);

ssim_map = filter2(window, img1, 'valid');        % gx
w1 = filter2(window, img2, 'valid');              % gy
w2 = ssim_map.*w1;                                % gx*gy
w2 = 2*w2+C1;                                     % 2*(gx*gy)+C1 = num1
w1 = (w1-ssim_map).^2+w2;                         % (gy-gx)^2+num1 = den1
ssim_map = filter2(window, img1.*img2, 'valid');  % g(x*y)
ssim_map = (2*ssim_map+(C1+C2))-w2;               % 2*g(x*y)+(C1+C2)-num1 = num2
ssim_map = ssim_map.*w2;                          % num
img1 = img1.^2;                                   % x^2
img2 = img2.^2;                                   % y^2
img1 = img1+img2;                                 % x^2+y^2

if (C1 > 0 && C2 > 0)
   w2 = filter2(window, img1, 'valid');           % g(x^2+y^2)
   w2 = w2-w1+(C1+C2);                            % den2
   w2 = w2.*w1;                                   % den
   ssim_map = ssim_map./w2;                       % num/den = ssim
else
   w3 = filter2(window, img1, 'valid');           % g(x^2+y^2)
   w3 = w3-w1+(C1+C2);                            % den2
   w4 = ones(size(w1));
   index = (w1.*w3 > 0);
   w4(index) = (ssim_map(index))./(w1(index).*w3(index));
   index = (w1 ~= 0) & (w3 == 0);
   w4(index) = w2(index)./w1(index);
   ssim_map = w4;
end

mssim = mean2(ssim_map);

return
