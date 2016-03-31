% REORIENT_FACETS Reorient faces of a triangle mesh (V,F) so that the left-hand
% rule normal of each face (consistently) points outward.
%
% [FF,I] = reorient_facets(V,F)
% [FF,I] = reorient_facets(V,F,'ParameterName',ParameterValue, ...)
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of triangle indicies into V
%   Optional:
%     'NumRays'  followed by total number of rays {#F*100}
%     'MinRays'  followed by minimum number of rays per patch/face {10}
%     'Facetwise'  followed by whether each facet should be considered
%       independently, could lead to inconsistent orientation of manifoldly
%       neighboring facets {false}
%     'UseParity'  Whether to use parity(?) {false}
% Outputs:
%   FF   #F by 3 list of reoriented facets
%   I  #F list of whether each face was flipped
%     
