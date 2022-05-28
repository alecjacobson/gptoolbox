function [C,alpha_C] = imdarken(A,alpha_A,B,alpha_B)
  % IMDARKEN Blend A over B, replacing pixels in B if the cooresponding pixel in
  % A is darker.
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

  % https://developer.android.com/reference/android/graphics/BlendMode#DARKEN
  alpha_C = alpha_A+alpha_B-alpha_B.*alpha_A;
  C = (1-alpha_B).*A + (1-alpha_A).*B+min(B,A);
end
