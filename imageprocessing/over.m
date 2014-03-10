function [C,alpha_C] = over(A,alpha_A,B,alpha_B)
  % OVER composite image A over image B using alpha masking
  %
  % [C,alpha_C] = over(A,alpha_A,B,alpha_B)
  %
  % Inputs:
  %   A  w by h by c image A
  %   alpha_A  alpha mask for image A, either:
  %     w by h by c
  %     w by h  (same across channels)
  %     s  (same for every pixel)
  %   B  w by h by c image B
  %   alpha_B  alpha mask for image B, either:
  %     w by h by c
  %     w by h  (same across channels)
  %     s  (same for every pixel)
  % Outputs:
  %   C  w by h by c result image, *ALWAYS of type double*
  %   alpha_C  alpha mask for resulting image
  %

  % only work with double images
  A = im2double(A);
  alpha_A = im2double(alpha_A);
  B = im2double(B);
  alpha_B = im2double(alpha_B);

  % composite images: http://en.wikipedia.org/wiki/Alpha_compositing
  C = bsxfun(@times,A,alpha_A) + bsxfun(@times,B,alpha_B.*(1-alpha_A));
  % composite alpha masks
  alpha_C = alpha_A + alpha_B.*(1-alpha_A);

end
