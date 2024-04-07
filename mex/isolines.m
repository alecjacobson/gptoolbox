% [iV,iE,I] = isolines(V,F,S,vals)
%
% Compute isolines of a scalar field on a triangle mesh.
%
% Isolines may cross perfectly at vertices. The output should not contain
% degenerate segments (so long as the input does not contain degenerate
% faces). The output segments are *oriented* so that isolines curl
% counter-clockwise around local maxima (i.e., for 2D scalar fields). Unless
% an isoline hits a boundary, it should be a closed loop. Isolines may run
% perfectly along boundaries. Isolines should appear just "above" constants
% regions.
% 
% Inputs:
%   V  #V by dim list of mesh vertex positions
%   F  #F by 3 list of mesh triangle indices into V
%   S  #S by 1 list of per-vertex scalar values
%   vals  #vals by 1 list of values to compute isolines for
% Outputs:
%   iV  #iV by dim list of isoline vertex positions
%   iE  #iE by 2 list of edge indices into iV
%   I  #iE by 1 list of indices into vals indicating which value each segment
%      belongs to

