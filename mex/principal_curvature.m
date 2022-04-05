% PRINCIPAL_CURVATURE Compute principal curvatures at vertices of a mesh (V,F) 
%
% [PD1,PD2,PV1,PV2] = principal_curvature(V,F);
%
% Inputs:
%   V  #V by 3 list of mesh vertex positions
%   F  #F by 3 list of mesh triangule indices into rows of V
% Outputs:
%   PD1  #V by 3 list of first principal curvature directions
%   PD2  #V by 3 list of second principal curvature directions
%   PV1  #V list of first principal curvature values
%   PV2  #V list of second principal curvature values
%
