% SIGNED_DISTANCE Compute signed distance from points P to a mesh (V,F)
%
% [S,I,C,N] = signed_distance(P,V,F,'ParameterName',parameter_value,...)
%
% Inputs:
%   P  #P by 3 list of query point positions
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of triangle indices
%   Optional:
%     'SignedDistanceType' followed by
%       'winding_number' use winding number (continuous sign value for
%         non-watertight)
%       'pseudonormal'  use pseudo-normal, binary scale (but not robust for
%         non-watertight meshes.
% Outputs:
%   S  #P list of smallest signed distances
%   I  #P list of facet indices corresponding to smallest distances
%   C  #P by 3 list of closest points
%   N  #P by 3 list of closest normals (only set if
%   
