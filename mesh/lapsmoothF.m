function [S] = lapsmoothF(TT,S,a,it)
% LAPSMOOTHF Laplacian smoothing for a scalar field defined on faces
% 
% [S] = lapsmoothF(TT,S,a,it)
% 
% Inputs:
%   TT  triangle-triangle topology (use tt to compute it)
%   S  scalar filed defined over faces, Fx1
%   a  controls the smoothing strenght (0..1)
%   it  number of iterations
% Output:
%   S  smoothed scalar field, Fx1

for x = 1:it
    for i=1:size(TT,1)
        acc = 0;
        count = 0;
        for j=1:3
            if(TT(i,j) ~= -1)
                acc = acc + S(TT(i,j));
                count = count + 1;
            end
        end
        S(i) = (1-a)*S(i) + ...
                (a * acc./count); 
    end
end
