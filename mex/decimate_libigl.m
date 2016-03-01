% DECIMATE_LIBIGL Decimate a closed manifold mesh (V,F)
%
% [W,G] = decimate_libigl(V,F,ratio)
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of triangle indices into V
%   ratio   either a 1<number<#F  of max faces, or a 0<ratio<1 to be multiplied
%     against #F to get max faces in output
% Outputs:
%   W  #W by 3 list of vertex positions
%   G  #G by 3 list of triangle indices into W
%
