function [VV] = vv(V,F)
% VV Compute the Vertex-Vertex topology of a manifold mesh
%
% [VV] = vv(V,F)
%
% Input:
%   V  vertex coordinates, Vx3
%   F  triangles         , Fx3
% Output:
%   VV  sparse matrix that contains a 1 in the i-j element if vertices i and j
%     are connected by edge

%% prepare sparse matrix
i = ones(10*size(V,1),1);
j = ones(10*size(V,1),1);
v = zeros(10*size(V,1),1);

row = 1;

for x=1:size(F,1) % for every face
    for y1=1:3 % for every pair of neighbours
        for y2=1:3
            if (y1 ~= y2)
                i(row) = F(x,y1);
                j(row) = F(x,y2);
                v(row) = 1;
                row = row + 1;
            end
        end
    end
end

VV = sparse(i,j,v);
end
