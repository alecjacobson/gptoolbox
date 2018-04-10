function [UV,F, res, edge_norms] = ...
  create_regular_grid(xRes, yRes, xWrap, yWrap, near, far)
% Creates list of triangle vertex indices for a rectangular domain,
% optionally wrapping around in X/Y direction.
%
% Usage:
%   [UV,F,res,edge_norms] = create_regular_grid(xRes, yRes, xWrap, yWrap)
%
% Input:
%   xRes, yRes: number of points in X/Y direction
%   wrapX, wrapY: wrap around in X/Y direction
%   near, far: near and far should be fractions of one which control the
%              pinching of the domain at the center and sides
%
% Output:
%   F : mesh connectivity (triangles)
%   UV: UV coordinates in interval [0,1]x[0,1]
%   res: mesh resolution
%
% Example:
%  % Create and m by n cylinder
%  m = 10; n = 20;
%  [V,F] = create_regular_grid(m,n,1,0);
%  V = [sin(2*pi*V(:,1)) cos(2*pi*V(:,1)) (n-1)*2*pi/(m-1)*V(:,2)];
%  tsurf(F,V); axis equal;
%  
%  % Quads:
%  Q = [F(1:2:end-1,[1 2]) F(2:2:end,[2 3])];
%

if (nargin<2) yRes=xRes; end
if (nargin<3) xWrap=0; end
if (nargin<4) yWrap=0; end
if (nargin<5) overlap=0; end

%res = [yRes, xRes];
res_wrap = [yRes+yWrap, xRes+xWrap];

%xSpace = linspace(0,1,xRes+xWrap); if (xWrap) xSpace = xSpace(1:end-1); end
%ySpace = linspace(0,1,yRes+yWrap); if (yWrap) ySpace = ySpace(1:end-1); end
xSpace = linspace(0,1,xRes+xWrap);
ySpace = linspace(0,1,yRes+yWrap);

[X, Y] = meshgrid(xSpace, ySpace);
UV_wrap = [X(:), Y(:)];

% Must perform pinch before edge_norms are taken
if(exist('near') & exist('far'))
  if(near>0 & far>0)
  t = ( ...
      UV_wrap(:,1).*(UV_wrap(:,1)<0.5)+ ...
      (1-UV_wrap(:,1)).*(UV_wrap(:,1)>=0.5) ...
    )/0.5;
  t = 1-sin(t*pi/2+pi/2);
  UV_wrap(:,2) = ...
    far/2 + ...
    near*(UV_wrap(:,2)-0.5).*(1-t) + ...
    far*(UV_wrap(:,2)-0.5).*t;
  else
    %error('Pinch must be between 0 and 1');
  end
end


idx_wrap = reshape(1:prod(res_wrap), res_wrap);

v1_wrap = idx_wrap(1:end-1, 1:end-1); v1_wrap=v1_wrap(:)';
v2_wrap = idx_wrap(1:end-1, 2:end  ); v2_wrap=v2_wrap(:)';
v3_wrap = idx_wrap(2:end  , 1:end-1); v3_wrap=v3_wrap(:)';
v4_wrap = idx_wrap(2:end  , 2:end  ); v4_wrap=v4_wrap(:)';

F_wrap = [v1_wrap;v2_wrap;v3_wrap; v2_wrap;v4_wrap;v3_wrap];
F_wrap = reshape(F_wrap, [3, 2*length(v1_wrap)])';

% old way
% edges = [F_wrap(:,1) F_wrap(:,2); F_wrap(:,2) F_wrap(:,3); F_wrap(:,3) F_wrap(:,1)];
% edge_norms = sqrt(sum((UV_wrap(edges(:,1),:)-UV_wrap(edges(:,2),:)).^2,2));
% edge_norms = reshape(edge_norms,size(F_wrap,1),3);

% edges numbered same as opposite vertices
edge_norms = [ ...
  sqrt(sum((UV_wrap(F_wrap(:,2),:)-UV_wrap(F_wrap(:,3),:)).^2,2)) ...
  sqrt(sum((UV_wrap(F_wrap(:,3),:)-UV_wrap(F_wrap(:,1),:)).^2,2)) ...
  sqrt(sum((UV_wrap(F_wrap(:,1),:)-UV_wrap(F_wrap(:,2),:)).^2,2)) ...
  ];

% correct indices
res = [yRes,xRes];
idx = reshape(1:prod(res),res);
if (xWrap) idx = [idx, idx(:,1)]; end
if (yWrap) idx = [idx; idx(1,:)]; end
idx_flat = idx(:);

% this might not be neccessary, could just rebuild UV like before
UV = reshape(UV_wrap,[size(idx_wrap),2]);
UV = UV(1:end-yWrap,1:end-xWrap,:);
UV = reshape(UV,xRes*yRes,2);

F = [idx_flat(F_wrap(:,1)),idx_flat(F_wrap(:,2)),idx_flat(F_wrap(:,3))];
