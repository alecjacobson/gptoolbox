function [G,GI] = fd_grad(side)
  % FD_GRAD build a finite difference gradient for a regular grid. Components
  % of the gradient live on **different** staggered grids. So, d/dx lives on a
  % grid staggered one half-step in the x-direction but otherwise aligned in
  % the y-direction (and z-direction).
  %
  % [G] = fd_grad(side)
  %
  % Inputs:
  %   dims  number of nodes along height (and width and depth)
  % Outputs:
  %   G   sum(prod(repmat(side,numel(side),1)-eye(numel(side)),2)) by
  %     prod(side)  gradient matrix [Gx;Gy; ...]
  %   GI   numel(side) array of row indices of Gx;Gy
  %
  % Example:
  %   G = fd_grad([3 3]);
  %   L = fd_laplacian([3 3]);
  %   -L - G'*G
  % 
  % Example:
  %   [X,Y] = meshgrid(linspace(-1,1,w),linspace(-1,1,h));
  %   X = X*(w-1)/(h-1);
  %   [G,GI] = fd_grad([w h]);
  %   Gx = G(GI(1)+1:GI(2),:);
  %   Gy = G(GI(2)+1:GI(3),:);
  %   % forward differences at bottom-left corners
  %   Gx(h:h:end,:) = [];
  %   Gy = Gy(1:end-(h-1),:);
  %   F = [Gx*Z(:) Gy*Z(:)];
  %   vec = @(X) X(:);
  %   Z = (X.^2+Y.^2);
  %   clf;
  %   hold on;
  %   surf(X,Y,0*Z,'CData',Z,fphong);
  %   quiver(vec(X(1:end-1,1:end-1)), vec(Y(1:end-1,1:end-1)),F(:,1),F(:,2),'k');
  %   hold off;
  %   axis equal;
  % 
  %
  % See also: fd_laplacian
  %

  function B = vec(A)
    B = A(:)';
  end


  G = [];
  n = prod(side);
  GI = [0];
  switch numel(side)
  case 2
    I = sub2ind([side(2) side(1)],1:n)';
    I = reshape(I,side([2 1]));
    for d = 1:2
      m = (side(2)-(d==2))*(side(1)-(d==1));
      J = sub2ind(  [side(2)-(d==2) side(1)-(d==1) ],1:m);
      J = reshape(J,[side(2)-(d==2) side(1)-(d==1) ]);
      G = [G;...
        sparse(vec(J),vec(I(1:end-(d==2),1:end-(d==1))),-1,m,n)+ ...
        sparse(vec(J),vec(I((d==2)+1:end,(d==1)+1:end)), 1,m,n)];
      GI = [GI size(G,1)];
    end 
  case 3
    I = sub2ind([side(2) side(1) side(3)],1:n)';
    I = reshape(I,side([2 1 3]));
    for d = 1:3
      m = (side(2)-(d==2))*(side(1)-(d==1))*(side(3)-(d==3));
      J = sub2ind(  [side(2)-(d==2) side(1)-(d==1) side(3)-(d==3)],1:m);
      J = reshape(J,[side(2)-(d==2) side(1)-(d==1) side(3)-(d==3)]);
      G = [G;...
        sparse(vec(J),vec(I(1:end-(d==2),1:end-(d==1),1:end-(d==3))),-1,m,n)+ ...
        sparse(vec(J),vec(I((d==2)+1:end,(d==1)+1:end,(d==3)+1:end)), 1,m,n)];
      GI = [GI size(G,1)];
    end 
  end 

end
