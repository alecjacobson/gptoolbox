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

  % expand (repmat) alphas to fit given images
  alpha_A = expand_alpha(alpha_A,size(A));
  alpha_B = expand_alpha(alpha_B,size(B));

  % composite images: http://en.wikipedia.org/wiki/Alpha_compositing
  C = A.*alpha_A + B.*alpha_B.*(1-alpha_A);
  % composite alpha masks
  alpha_C = alpha_A + alpha_B.*(1-alpha_A);

  function alpha = expand_alpha(alpha,S)
    % EXPAND_ALPHA expand alpha mask to fit given size by repmating for
    % channels or width,height,channels
    % 
    % alpha = expand_alpha(alpha,S)
    % 
    % Inputs
    %   alpha  given alpha mask
    %   S  size of image corresponding to alpha mask
    % Outputs
    %   alpha  alpha map resized so that size(alpha) == S
    %
    if numel(size(alpha)) == numel(S) && all(size(alpha) == S)
      % do nothing
    elseif numel(alpha) == 1
      % repmat for all pixels
      alpha = repmat(alpha,S);
    elseif size(alpha,1) == S(1) && size(alpha,2) == S(2) && numel(S) == 3
      % just repmat for channels
      alpha = repmat(alpha,[1 1 S(3)]);
    else
      error('Given alpha mask is not compatible with given size');
    end
  end
end
