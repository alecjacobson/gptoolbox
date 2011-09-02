function [UV,F, res] = create_irregular_grid(xRes, yRes, n, xWrap, yWrap)
% creates random UV positions and their connectivity information.
%
% Usage:
%   [UV,F,res] = create_irregular_grid(xRes, yRes, n, wrapX, wrapY)
%
% Input:
%    n: number of points to create per block
%    wrapX, wrapY: wrap around in X/Y direction
% Output:
%   UV: UV coordinates in interval [0,1]x[0,1]
%   F : mesh connectivity (triangles)

res = [yRes, xRes];

xSpace = linspace(0,1,xRes+1); xSpace=xSpace(1:end-1);
ySpace = linspace(0,1,yRes+1); ySpace=ySpace(1:end-1);
uvR = rand(xRes*yRes*n, 2).*repmat([xSpace(2), ySpace(2)], [xRes*yRes*n, 1]);
[U,V] = meshgrid(xSpace, ySpace);
uvB = repmat([U(:), V(:)], [n 1]);

uvB(uvB(:,1)<0.001,1) = 0.001;
uvB(uvB(:,2)<0.001,2) = 0.001;
uvB(uvB(:,1)>0.999,1) = 0.999;
uvB(uvB(:,2)>0.999,2) = 0.999;

uv = uvR + uvB;
nUV = size(uv,1);
nX = round(sqrt(n)*xRes); nY = round(sqrt(n)*xRes);
xBorder = [0; rand(nX, 1); 1];
yBorder = [rand(nY, 1)];

UV = [uv; [xBorder, 0*xBorder]; [0*yBorder, yBorder];
    [xBorder, 0*xBorder+1]; [0*yBorder+1, yBorder]];

nUV = nUV + nX + nY;

F = delaunay(UV(:,1), UV(:,2));

if (xWrap)
    if (yWrap)
        F(F>nUV) = F(F>nUV) - nX - nY;
        UV = UV(1:nUV, :);
    else
        F(F>nUV & F<=nUV+nX) = F(F>nUV & F<=nUV+nX) - nX - nY;
        F(F>nUV) = F(F>nUV) - nX;
        UV = [UV(1:nUV, :), [0*yBorder+1, yBorder]]; 
    end
%else
%    nUV = nUV + nX;
%    F(F>nUV) = F(F>nUV) - nX - nY;
%    UV = UV(1:nUV,:);
end




end
