function str = falpha(alpha,ealpha)
  % FPHONG  Simply returns:
  % struct('FaceAlpha',alpha,'EdgeAlpha',alpha);
  %
  % str = fphong()
  % 
  if nargin == 0
    alpha = 0.6;
  end
  if nargin == 1
    ealpha = alpha;
  end
  str = struct('FaceAlpha',alpha,'EdgeAlpha',ealpha);
end

