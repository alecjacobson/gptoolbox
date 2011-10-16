function [ SV ] = moveFV(V,F,S)
% moveFV 
% Move a scalar field defined on faces to vertices by averaging
%
% Input:
% V,F: mesh
% S: scalar field defined on faces, Fx1
% 
% Output:
% SV: scalar field defined on vertices

SV = zeros(size(V,1),1);
COUNT = zeros(size(V,1),1);

for i=1:size(F,1)
    SV(F(i,:)') = SV(F(i,:)') + S(i);
    COUNT(F(i,:)') = COUNT(F(i,:)') + 1;
end

SV = SV ./ COUNT;

end

