function expand_axis(f,ah)
  % EXPAND_AXIS Expand axis by a scale factor f
  % 
  % expand_axis(f);
  % expand_axis(f,ah)
  %
  % Inputs:
  %   f  factor to scale axis by
  %   Optional:
  %   ah  handle to axis
  %  
  if nargin<2
    ah =gca;
  end
  a = reshape(axis(ah),2,[]);
  axis(ah,reshape(bsxfun(@plus,bsxfun(@minus,a,mean(a))*f,mean(a)),1,[]));
end
