function Z = rasterize(V,E)
  % RASTERIZE A weirdly vectorized implementation of Bresenham-like
  % rasterization of 2D polylines
  %
  % Z = rasterize(V,E)
  %
  % Inputs:
  %   V  #V by 2 list of vertex positions
  %   E  #E by 2 list of edge indices into rows of V
  % Outputs:
  %   w by h  "image" of indices into rows of E (0 means pixel is not drawn
  %     upon)
  %


  % expect that coordinates are already positive integers
  x0 = round(V(E(:,1),1));
  y0 = round(V(E(:,1),2));
  x1 = round(V(E(:,2),1));
  y1 = round(V(E(:,2),2));

  % Z [1 1] to ceil(max(V))
  assert(floor(min(V(:)))>=1);
  sqrlen = ((x0-x1).^2 + (y0-y1).^2);
  II = find(sqrlen>0);
  x0 = x0(II);
  x1 = x1(II);
  y0 = y0(II);
  y1 = y1(II);

  Z = zeros(fliplr(ceil(max(V))));

  % x increasing
  [Xx,I] = sort([x0 x1],2);
  Xy = [y0 y1];
  Xy(I(:,1)==2,[1 2]) = Xy(I(:,1)==2,[2 1]);
  Xlen = Xx(:,2)-Xx(:,1)+1;

  % y increasing
  [Yy,I] = sort([y0 y1],2);
  Yx = [x0 x1];
  Yx(I(:,1)==2,[1 2]) = Yx(I(:,1)==2,[2 1]);
  Ylen = Yy(:,2)-Yy(:,1)+1;

  % More horizontal lines
  I = Xlen>Ylen;
  Ix = Xx(I,:);
  Iy = Xy(I,:);
  Ilen = Xlen(I);
  %Ix = repelem(Ix,Ilen)
  IX = cumsum(repelem(ones(numel(Ilen),1),Ilen)) - ...
    repelem(cumsum([0;Ilen(1:end-1)]),Ilen)-1+ ...
    repelem(Ix(:,1),Ilen);
  Id = (Iy(:,2)-Iy(:,1))./(Ix(:,2)-Ix(:,1));
  IY = cumsum(repelem(Id,Ilen))-repelem(cumsum([0;Id(1:end-1).*(Ilen(1:end-1))]),Ilen)+repelem(Iy(:,1),Ilen) - repelem(Id,Ilen);
  IY = round(IY);
  Z(sub2ind(size(Z),IY,IX)) = repelem(find(I),Ilen);

  % More vertical lines
  I = Xlen<=Ylen;
  Iy = Yy(I,:);
  Ix = Yx(I,:);
  Ilen = Ylen(I);
  %Iy = repelem(Iy,Ilen)
  IY = cumsum(repelem(ones(numel(Ilen),1),Ilen)) - ...
    repelem(cumsum([0;Ilen(1:end-1)]),Ilen)-1+ ...
    repelem(Iy(:,1),Ilen);
  Id = (Ix(:,2)-Ix(:,1))./(Iy(:,2)-Iy(:,1));
  IX = cumsum(repelem(Id,Ilen))-repelem(cumsum([0;Id(1:end-1).*(Ilen(1:end-1))]),Ilen)+repelem(Ix(:,1),Ilen) - repelem(Id,Ilen);
  IX = round(IX);

  Z(sub2ind(size(Z),IY,IX)) = repelem(II(find(I)),Ilen);


  %IX = repelem(Id,Ilen) + ...
  %  repelem(cumsum([0;Id(1:end-1).*Ilen(1:end-1)]),Ilen)-1+ ...




end
