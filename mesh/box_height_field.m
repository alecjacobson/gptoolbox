function [V,F] = box_height_field(im)
  % Input:
  %   im  h by w by 1 height field over regular grid
  % Outputs:
  %   V  #V by 3 list of 3d surface mesh positions
  %   F  #F by 3 list of 3d surface triangle indices
  %

  function [V,F] = strips(X,Y,Z)
    %
    %         TL2-TR2-...
    %        / |   |
    % TL1-TR1 BL2-BR2-...
    %  |   | /
    % BL1+BR1
    %
    BLX = X-0.5;
    BRX = X+0.5;
    BLY = Y-0.5;
    BRY = Y-0.5;
    BLZ = Z;
    BRZ = Z;
    TLX = X-0.5;
    TRX = X+0.5;
    TLY = Y+0.5;
    TRY = Y+0.5;
    TLZ = Z;
    TRZ = Z;
    BLI = reshape(0*numel(X)+(1:numel(X)),size(X));
    BRI = reshape(1*numel(X)+(1:numel(X)),size(X));
    TLI = reshape(2*numel(X)+(1:numel(X)),size(X));
    TRI = reshape(3*numel(X)+(1:numel(X)),size(X));
    V = [ ...
      BLX(:) BLY(:) BLZ(:); ...
      BRX(:) BRY(:) BRZ(:); ...
      TLX(:) TLY(:) TLZ(:); ...
      TRX(:) TRY(:) TRZ(:)];
    F = [ ...
      BLI(:) BRI(:) TRI(:);
      BLI(:) TRI(:) TLI(:);
      reshape(BRI(:,1:end-1),[],1), ...
        reshape(BLI(:,2:end),[],1), ...
        reshape(TLI(:,2:end),[],1); ...
      reshape(BRI(:,1:end-1),[],1), ...
        reshape(TLI(:,2:end),[],1), ...
        reshape(TRI(:,1:end-1),[],1); ...
      ];
  end

  % Gather X,Y,Z for each pixel
  w = size(im,2);
  h = size(im,1);
  [X,Y] = meshgrid(1:w,1:h);
  Z = im;
  %X = X(:);
  %Y = Y(:);
  %Z = im(:);
  I = reshape(1:numel(im),h,w);

  % left-right strips
  [VLR,FLR] = strips(X,Y,Z);
  [VTB,FTB] = strips(Y',X',Z');
  VTB = VTB(:,[2 1 3]);
  FTB = fliplr(FTB);
  V = [VLR;VTB];
  F = [FLR;size(VLR,1)+FTB];

  [V,F] = clean(V,F,'MinDist',eps,'MinArea',0,'MinAngle',0, ...
       'SelfIntersections','mesh','SmallTriangles','remove');
  % This can result in a non-manifold mesh if there is every a neighborhood of
  % 4 pixels creating a "saddle" at a corner. For example,
  %
  %  0 1
  %  1 0
  %
  %  Then the edge running _up_ the corner between the two 1 towers will be
  %  non-manifold.
  %  

  %%% Create 8 vertices for each pixel
  %%%
  %%%        TD--TC
  %%%       /|    |
  %%%      / TA--TB
  %%%     /  /    /
  %%%    /  /    /
  %%%   /  /    /
  %%%  BD--BC /
  %%%  |/   |/
  %%%  BA--BB
  %%%
  %%%

  %%BAX = X(2:end,2:end)-0.5;
  %%BAY = Y(2:end,2:end)-0.5;
  %%BAZ = Z(1:end-1,1:end-1);

  %%BBX = X(2:end,2:end)-0.5;
  %%BBY = Y(2:end,2:end)-0.5;
  %%BBZ = Z(1:end,1:end);

end
