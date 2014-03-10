function [ SV ] = moveFV(V,F,S)
% MOVEFV  Move a scalar field defined on faces to vertices by averaging
% 
% [ SV ] = moveFV(V,F,S)
%
% Input:
%   V,F  mesh
%   S  scalar field defined on faces, Fx1
% 
% Output:
%   SV  scalar field defined on vertices

SV = zeros(size(V,1),size(S,2));
COUNT = zeros(size(V,1),1);

for i=1:size(F,1)
    SV(F(i,:)',:) = SV(F(i,:)',:) + repmat(S(i,:),size(F,2),1);
    COUNT(F(i,:)') = COUNT(F(i,:)') + 1;
end

SV = SV ./ repmat(COUNT, 1, size(S,2));

end

