% DUAL_LAPLACIAN from "Laplace Operators for Tetrahedral Meshes" [Alexa et al.
% 2020]
% 
% [L,M] = dual_laplacian(V,T)
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   T  #T by 4 list of tetrahedron indices into rows of V
% Outputs:
%   L  #V by #V sparse dual Laplacian matrix
%   M  #V by #V sparse mass matrix
