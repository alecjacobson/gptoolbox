% SIGNED_DISTANCE_ISOSURFACE Compute the contour of an iso-level of the
% signed distance field to a given mesh.
%
% [V,F] = signed_distance_isosurface(IV,IF);
% [V,F] = signed_distance_isosurface(IV,IF,'ParameterName',parameter_value, ...);
%
% Inputs:
%   IV  #IV by 3 list of input mesh vertex positions
%   IF  #IF by 3 list of input triangle indices
%   Optional:
%     'Level' followed by iso-level to contour in world coordinate units {0}
%     'AngleBound' followed by lower bound on triangle angles (mesh quality)
%       {28}
%     'RadiusBound' followed by upper bound on triangle size (mesh density?)
%       as fraction of bounding box diagonal {0.02} 
%     'DistanceBound' followed by cgal mysterious parameter (mesh density?)
%       as fraction of bounding box diagonal {0.02} 
% Outputs:
%   V  #V by 3 list of output mesh positions
%   F  #F by 3 list of output triangle indices into V
%

