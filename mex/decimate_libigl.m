% DECIMATE_LIBIGL Decimate a closed manifold mesh (V,F)
%
% [W,G] = decimate_libigl(V,F,ratio)
% [W,G,J,I] = decimate_libigl(V,F,ratio,'ParameterName',ParameterValue, ...)
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of triangle indices into V
%   ratio   either a 1<number<#F  of max faces, or a 0<ratio<1 to be multiplied
%     against #F to get max faces in output
%   Optional:
%     'Method' followed by one of:
%        {'naive'}  simply collapse small edges and place vertices at midpoint
%        'qslim'  Quadric error metric
% Outputs:
%   W  #W by 3 list of vertex positions
%   G  #G by 3 list of triangle indices into W
%   J  #G list of indices into F of birth face
%   I  #U list of indices into V of birth vertices
%
