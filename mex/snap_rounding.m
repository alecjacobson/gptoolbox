% SNAP_ROUNDING Snap a list of possible intersecting segments with
% endpoints in any precision to _the_ integer grid.
%
% Inputs:
%   V  #V by 2 list of vertex positions
%   E  #E by 2 list of segment indices into V
% Outputs:
%   VI  #VI by 2 list of output integer vertex positions, rounded copies
%     of V are always the first #V vertices
%   EI  #EI by 2 list of segment indices into V, #EI ≥ #E
%   J  #EI list of indices into E revealing "parent segments"
%
