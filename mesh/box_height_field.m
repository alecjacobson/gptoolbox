function [V,F] = box_height_field(im)
  % Create a height field of boxes (each pixel is an extruded square). The
  % resulting mesh will be water-tight (the surface of the solid beneath the
  % height field) except for the outer boundary. The surface will likely _not_
  % be edge manifold (saddle corners with data [1 0;0 1] will result in
  % non-manifold edges along the kissing pillars).
  %
  % Input:
  %   im  h by w by 1 height field over regular grid
  % Outputs:
  %   V  #V by 3 list of 3d surface mesh positions
  %   F  #F by 3 list of 3d surface triangle indices
  %

  function [V,F] = strips(X,Y,Z)
    F = [];
    V = [];
    h = size(X,1);
    w = size(X,2);
    for r = 1:h
      %              TL2-TR2-...
      %             / |   |
      %            / BL2-BR2-...
      %           / /
      %          / /
      %         /|/ 
      %        / ^ 
      % TL1-TR1 /
      %  |   | /
      % BL1+BR1
      %
      %% get next and previous row (or same if boundary)
      %Z2 = im(min(r+1,end),1);
      %Z1 = im(r,1);
      %Z0 = im(max(r-1,1),1);
      Vrc = [ ...
        X(r,1)-0.5 Y(r,1)-0.5 Z(r,1); ...
        X(r,1)-0.5 Y(r,1)+0.5 Z(r,1); ...
        X(r,1)+0.5 Y(r,1)+0.5 Z(r,1); ...
        X(r,1)+0.5 Y(r,1)-0.5 Z(r,1); ...
        ];
      Frc = [1 2 4;4 2 3];
      F = [F;size(V,1)+Frc];
      V = [V;Vrc];
      for c = 1:w-1
        r_prev = max(r,1);
        r_next = min(r,h);
        Zabove = zeros(0,1);
        Zbelow = zeros(0,1);
        if r>1
          Zabove = Z(r-1,[c c+1])';
          Zabove = sort(Zabove( ...
            (Zabove>Z(r,c) & Zabove<Z(r,c+1)) | ...
            (Zabove<Z(r,c) & Zabove>Z(r,c+1))),'descend');
        end
        if r<h
          Zbelow = Z(r+1,[c c+1])';
          Zbelow = sort(Zbelow( ...
            (Zbelow>Z(r,c) & Zbelow<Z(r,c+1)) | ...
            (Zbelow<Z(r,c) & Zbelow>Z(r,c+1))));
        end
        if Z(r,c) > Z(r,c+1)
          Zbelow = flipud(Zbelow);
          Zabove = flipud(Zabove);
        end
    
        Vrc = [ ...
          X(r,c)+0.5 Y(r,c)-0.5 Z(r,c); ...
          X(r,c)+0.5 Y(r,c)+0.5 Z(r,c); ...
          repmat([X(r,c)+0.5 Y(r,c)+0.5],numel(Zbelow),1) Zbelow; ...
          X(r,c)+0.5 Y(r,c)+0.5 Z(r,c+1); ...
          X(r,c)+0.5 Y(r,c)-0.5 Z(r,c+1); ...
          repmat([X(r,c)+0.5 Y(r,c)-0.5],numel(Zabove),1) Zabove; ...
          ];
        nbelow = 2+numel(Zbelow)+1;
        nabove = 2+numel(Zabove)+1;
        Frc = [ ...
          [ones(1,nbelow-2); 2:nbelow-1; 3:nbelow]'; ...
          nbelow-1+[ones(1,nabove-2); 2:nabove-1; [3:nabove-1 -nbelow+2]]'];
    
    
        Frc = [Frc; ...
          [1 [2 1]-2-numel(Zabove); ...
          2 1 [1]-2-numel(Zabove); ...
          ]+size(Vrc,1)];
        Vrc = [Vrc; ...
          X(r,c+1)+0.5 Y(r,c+1)-0.5 Z(r,c+1); ...
          X(r,c+1)+0.5 Y(r,c+1)+0.5 Z(r,c+1)];
    
        F = [F;size(V,1)+Frc];
        V = [V;Vrc];
      end
    end
  end

  w = size(im,2);
  h = size(im,1);
  [X,Y] = meshgrid(1:w,1:h);
  Z = im;
  [VLR,FLR] = strips(X,Y,Z);
  [VTB,FTB] = strips(Y',X',Z');
  VTB = VTB(:,[2 1 3]);
  V = [VLR;VTB];
  F = [FLR;size(VLR,1)+fliplr(FTB)];

  % clean up
  [V,~,IM] = remove_duplicate_vertices(V,0);
  F = IM(F);
  [~,I] = unique(sort(F,2),'rows');
  F = F(I,:);
  F = F(doublearea(V,F)>0,:);
  [V,IM] = remove_unreferenced(V,F);
  F = IM(F);

end
