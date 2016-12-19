function A = uninterp(V,BC,side)
  % UNINTERP un-interpolate some scattered data at points V on a grid defined
  % by (BC,side). That is construct an interpolation matrix A, so that each row
  % A(i,:) corresponds to an interpolation equation:
  %
  %  A(i,:) * BCZ(:) = Z(i)
  %
  % where BCZ would be values stored on the (BC,side) grid and Z(i) would be
  % the value at the ith data point at location V(i,:)
  %
  % Inputs:
  %   V  #V by dim list of query point locations
  %   BC  prod(side) by dim list of grid node positions
  %   side  dim-long list of nodes in x-direction, y-direction, ...
  % Outputs:
  %   A  #q by size(BC,1) sparse matrix
  %
  % Example:
  %
  % See also: interp2, interp3, interpn, meshgrid
  %

  dim = size(V,2);
  assert(numel(side) == dim);
  assert(size(BC,1) == prod(side));

  switch dim
  case 2
    Xq = V(:,1);
    Yq = V(:,2);
    X = reshape(BC(:,1),side([2 1]));
    Y = reshape(BC(:,2),side([2 1]));
    % Xq  #q by 1 list of x-coordinates of scatter data 
    % Yq  #q by 1 list of y-coordinates of scatter data 
    % X  x-coordinates of grid (like meshgrid)
    % Y  y-coordinates of grid (like meshgrid)
    [~,x] = histc(Xq,X(1,:));
    x(x==size(X,2)) = x(x==size(X,2))-1;
    [~,y] = histc(Yq,Y(:,1));
    y(y==size(Y,1)) = y(y==size(Y,1))-1;
    tx = (Xq-X(1,x)')./(X(1,x+1)'-X(1,x)');
    ty = (Yq-Y(y,1) )./(Y(y+1,1) -Y(y,1) );
    A = sparse( ...
     repmat(1:size(Xq,1),4,1)', ...
     sub2ind(size(X),[y y y+1 y+1],[x x+1 x+1 x]), ...
     [(1-tx).*(1-ty) tx.*(1-ty) tx.*ty (1-tx).*ty], ...
     size(Xq,1),numel(X));
  case 3
    Xq = V(:,1);
    Yq = V(:,2);
    Zq = V(:,3);
    X = reshape(BC(:,1),side([2 1 3]));
    Y = reshape(BC(:,2),side([2 1 3]));
    Z = reshape(BC(:,3),side([2 1 3]));
    [~,x] = histc(Xq,X(1,:,1));
    x(x==size(X,2)) = x(x==size(X,2))-1;
    [~,y] = histc(Yq,Y(:,1,1));
    y(y==size(Y,1)) = y(y==size(Y,1))-1;
    [~,z] = histc(Zq,Z(1,1,:));
    z(z==size(Z,3)) = z(z==size(Z,3))-1;
    tx = squeeze((permute(Xq,[2 1 3])-X(1,x,1))./(X(1,x+1,1)-X(1,x,1)))';
    ty = squeeze((permute(Yq,[1 2 3])-Y(y,1,1))./(Y(y+1,1,1)-Y(y,1,1)));
    tz = squeeze((permute(Zq,[3 2 1])-Z(1,1,z))./(Z(1,1,z+1)-Z(1,1,z)));
    A = sparse( ...
     repmat(1:size(Xq,1),8,1)', ...
     sub2ind(size(X), ...
       [y   y   y+1 y+1 y   y   y+1 y+1], ...
       [x   x+1 x+1 x   x   x+1 x+1 x  ], ...
       [z   z   z   z   z+1 z+1 z+1 z+1]), ...
       [(1-tx).*(1-ty).*(1-tz) ...
            tx.*(1-ty).*(1-tz) ...
                tx.*ty.*(1-tz) ...
                    (1-tx).*ty.*(1-tz) ...
                        (1-tx).*(1-ty).*tz ...
                            tx.*(1-ty).*tz ...
                                tx.*ty.*tz ...
                                    (1-tx).*ty.*tz], ...
     size(Xq,1),numel(X));
  otherwise
    error('Dimension not supported: %d',dim);
  end

end
