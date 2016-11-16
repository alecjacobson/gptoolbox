function [G] = fd_grad(side)
  % FD_GRAD build a finite difference gradient for a regular grid. Components
  % of the gradient live on **different** staggered grids. So, d/dx lives on a
  % grid staggered one half-step in the x-direction but otherwise aligned in
  % the y-direction (and z-direction).
  %
  % Inputs:
  %   side  containing number of vertices on each side of grid [x y z]
  % Outputs:
  %   G   sum(prod(repmat(side,numel(side),1)-eye(numel(side)),2)) by
  %     prod(side)  gradient matrix [Gx;Gy; ...]
  %
  % Example:
  %   G = fd_grad([3 3]);
  %   L = fd_laplacian([3 3]);
  %   -L - G'*G
  %
  % See also: fd_laplacian
  %

  function B = vec(A)
    B = A(:)';
  end


  G = [];
  n = prod(side);
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
    end 
  end 

end
