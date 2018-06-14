function [A,C] = fd_trilinear_coefficients(mn,mx,side,V)
  % FD_TRILINEAR_COEFFICIENTS  Given a grid from minimum corner mn to maximum
  % corner mx, with `side` nodes on each side, compute a matrix `A` so that
  % `V=A*C` if `C` are the node centers.
  %
  % [A,C] = fd_trilinear_coefficients(mn,mx,side,V)
  %
  % Inputs:
  %   mn  3d position of minimum corner
  %   mx  3d position of maximum corner
  %   side  number of nodes in each dimension [x y z]
  %   V  #V by 3 list of query points
  % Outputs:
  %   A  x*y*z by #V matrix
  %   C  x*y*z by 3 list of node center positions 
  %   

  [X,Y,Z] = meshgrid( ...
    linspace(mn(1),mx(1),side(1)), ...
    linspace(mn(2),mx(2),side(2)), ...
    linspace(mn(3),mx(3),side(3)));
  C = [X(:) Y(:) Z(:)];
  r = (mx-mn)./(side-1);

  I = floor(bsxfun(@times,bsxfun(@rdivide,bsxfun(@minus,V,mn),mx-mn),side-1));
  A = sparse(size(V,1),prod(side));
  for x = 1:2
    for y = 1:2
      for z = 1:2
        Ixyz = sub2ind(side([2 1 3]),I(:,2)+x,I(:,1)+y,I(:,3)+z);
        Sxyz = C(Ixyz,:);
        Axyz = prod(bsxfun(@minus,r,abs(V-Sxyz)),2)/prod(r);
        A = A + sparse((1:size(V,1))',Ixyz,Axyz,size(V,1),size(C,1));
      end
    end
  end

end
