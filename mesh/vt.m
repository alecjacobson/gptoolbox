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

%% prepare sparse matrix
i = ones(10*size(V,1),1);
j = ones(10*size(V,1),1);
v = zeros(10*size(V,1),1);

row = 1;

for x=1:size(F,1) % for every face
    for y1=1:3 % for every pair of neighbours
        i(row) = x;
        j(row) = F(x,y1);
        v(row) = 1;
        row = row + 1;
    end
end

VT = sparse(i,j,v);
end
