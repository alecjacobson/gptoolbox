% FIT_CUBIC_BEZIER Fit a cubic bezier spline (G1 continuous) to an ordered list
% of input points in any dimension, according to "An algorithm for automatically
% fitting digitized curves" [Schneider 1990].
% 
%  Inputs:
%    d  #d by dim list of points along a curve to be fit with a cubic bezier
%      spline (should probably be roughly uniformly spaced). If d(0)==d(end),
%      then will treat as a closed curve.
%    error  maximum squared distance allowed
%  Output:
%    cubics #cubics list of 4 by dim lists of cubic control points
%
