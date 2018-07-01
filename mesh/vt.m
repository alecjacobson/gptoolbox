function [VT] = vt(V,F)
% VT Compute the Vertex-Face topology
%
% VT = vt(V,F)
%
% Input:
%   V  vertex coordinates, Vx3
%   F  triangles         , Fx3
% Output:
%   VT  sparse matrix that contains a 1 in the i-j element if vertices i is
%     part of face j

i = (1:size(F,1))';
j = F;
VT = sparse([i i i],j,1);
end
