% MESH_BOOLEAN Compute boolean csg operations on "solid", consistently oriented
% meshes.
%
% [W,H] = mesh_boolean(V,F,U,G,operation)
% [W,H] = mesh_boolean(V,F,U,G,operation,'ParameterName',paramter_value, ...)
% 
% Inputs:
%   V  #V by 3 list of vertex positions of first mesh
%   F  #F by 3 list of triangle indices into V
%   U  #U by 3 list of vertex positions of second mesh
%   G  #G by 3 list of triangle indices into U
%   operation  followed by operation to perform as a string, one of: 'union',
%     'intersect', 'minus', 'xor', or 'resolve'
%     Optional:
%       'BooleanLib' followed by boolean library back-end to use, one of:
%         {'libigl'}  uses CGAL's exact arithmetic kernel and is believed to be
%                     correct.
%         'cork'  is faster but may give incorrect results. 
%         'libigl-try-cork-resolve'  libigl boolean extraction but tries to use
%                                    cork's fast resolve, if intersections
%                                    persist, then resolves remaining with
%                                    libigl's resolve. This adds a "layer of
%                                    robustness" on top of cork, but since it's
%                                    not understood _how_ cork is failing, it
%                                    is unknown whether this will lead to
%                                    correct results.
% Outputs:
%   W  #W by 3 list of vertex positions of boolean result mesh
%   H  #H by 3 list of triangle indices into W
%   J  #H list of indices into [FA;FB] of facet birth parents
% 
% See also: self_intersect
%    
