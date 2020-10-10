function [A,C] = fd_bilinear_coefficients(mn,mx,side,V)
  % FD_BILINEAR_COEFFICIENTS  Given a grid from minimum corner mn to maximum
  % corner mx, with `side` nodes on each side, compute a matrix `A` so that
  % `V=A*C` if `C` are the node centers.
  %
  % [A,C] = fd_bilinear_coefficients(mn,mx,side,V)
  %
  % Inputs:
  %   mn  3d position of minimum corner
  %   mx  3d position of maximum corner
  %   side  number of nodes in each dimension [x y]
  %   V  #V by 2 list of query points
  % Outputs:
  %   A  x*y by #V matrix
  %   C  x*y by 3 list of node center positions 
  %   

  [X,Y] = meshgrid( ...
    linspace(mn(1),mx(1),side(1)), ...
    linspace(mn(2),mx(2),side(2)));
  C = [X(:) Y(:)];
  r = (mx-mn)./(side-1);

  I = floor(bsxfun(@times,bsxfun(@rdivide,bsxfun(@minus,V,mn),mx-mn),side-1));
  A = sparse(size(V,1),prod(side));
  for x = 1:2
    for y = 1:2
      Ixy = sub2ind(side([2 1]),I(:,2)+x,I(:,1)+y);
      Sxy = C(Ixy,:);
      Axy = prod(bsxfun(@minus,r,abs(V-Sxy)),2)/prod(r);
      A = A + sparse((1:size(V,1))',Ixy,Axy,size(V,1),size(C,1));
    end
  end

end
