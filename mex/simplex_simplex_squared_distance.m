% SIMPLEX_SIMPLEX_SQUARED_DISTANCE Squared distance between two simplices.
%
% [sqrd,C1,C2] = simplex_simplex_squared_distance(V1,V2)
% [sqrd,C1,C2,B1,B2] = simplex_simplex_squared_distance(V1,V2)
%
% Inputs:
%   V1  #V1 by dim list of corners of a (#V1-1)-simplex
%   V2  #V2 by dim list of corners of a (#V2-1)-simplex
% Outputs:
%   sqrd  squared distance between the simplices
%   C1  1 by dim closest point on the first simplex
%   C2  1 by dim closest point on the second simplex
%   B1  1 by #V1 barycentric coordinates of C1 w.r.t. V1 (so C1 == B1*V1)
%   B2  1 by #V2 barycentric coordinates of C2 w.r.t. V2 (so C2 == B2*V2)
%
% The simplices may have different numbers of corners but must have the same
% dimension. Degenerate simplices are handled.
%
% Example:
%   % squared distance between two triangles in 3D
%   [sqrd,C1,C2] = simplex_simplex_squared_distance(V1(F1(i,:),:),V2(F2(j,:),:));
%
% See also: point_simplex_squared_distance, aabb_aabb_squared_distance,
% mesh_mesh_squared_distance
%
