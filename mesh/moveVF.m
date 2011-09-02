function [ SF ] = moveVF(F,S)
% moveVF 
% Move a scalar field defined on vertices to faces by averaging
%
% Input:
% V,F: mesh
% S: scalar field defined on vertices, Vx1
% 
% Output:
% SF: scalar field defined on faces

M = sparse(repmat((1:size(F,1))',3,1),[F(:,1),F(:,2),F(:,3)],1/3);
SF = M*S;

end

