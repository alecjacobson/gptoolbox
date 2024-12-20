% EXACT_GEODESIC
%
% Inputs:
%   V  #V by dim list of vertex positions
%   F  #F by 3 list of triangle indices into rows of V
%   VS  #VS by 1 list of indices into rows of V of source vertices
%   FS  #FS by 1 list of indices into rows of F of source faces
%   VT  #VT by 1 list of indices into rows of V of target vertices
%   FT  #FT by 1 list of indices into rows of F of target faces
% Outputs:
%   D  #VT+#FT by 1 list of geodesic distances of each target w.r.t. the nearest one in the source set
%
%
