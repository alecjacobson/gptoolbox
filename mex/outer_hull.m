% OUTER_HULL Compute the "outer hull" of a potentially non-manifold mesh (V,F)
% whose intersections have been "resolved" (e.g. using `cork` or
% `igl::selfintersect`). The outer hull is defined to be all facets (regardless
% of orientation) for which there exists some path from infinity to the face
% without intersecting any other facets. For solids, this is the surface of the
% solid. In general this includes any thin "wings" or "flaps".  This
% implementation largely follows Section 3.6 of "Direct repair of
% self-intersecting meshes" [Attene 2014].
%
% [G,J,flip] = outer_hull(V,F);
%
% Inputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of triangle indices into V
% Outputs:
%   G  #G by 3 list of output triangle indices into V
%   J  #G list of indices into F
%   flip  #F list of whether facet was added to G **and** flipped orientation
%     (false for faces not added to G)
