function [V,F] = box_height_field(im,varargin)
  % Create a height field of boxes (each pixel is an extruded square). The
  % resulting mesh will be water-tight (the surface of the solid beneath the
  % height field) except for the outer boundary. The surface will likely _not_
  % be edge manifold (saddle corners with data [1 0;0 1] will result in
  % non-manifold edges along the kissing pillars).
  %
  % [V,F] = box_height_field(im)
  % [V,F] = box_height_field(im,'ParameterName',ParameterValue, ...)
  %
  % Input:
  %   im  h by w by 1 height field over regular grid
  %   Optional:
  %     'Closed' whether to produce a closed mesh (slower) {true}.
  % Outputs:
  %   V  #V by 3 list of 3d surface mesh positions
  %   F  #F by 3 list of 3d surface triangle indices
  %

  function [V,F] = box_height_field_closed(im)
    function [V,F] = strips(X,Y,Z,tops)
      h = size(X,1);
      w = size(X,2);
      % reserve space
      F = zeros(h*w*6,3);
      V = zeros(h*w*10,3);
      fsize = 0;
      vsize = 0;
  
      function appendF(Frc)
        if fsize+size(Frc,1) > size(F,1)
          Ftemp = F;
          F = [F;zeros(size(F,1),3)];
        end
        F(fsize+(1:size(Frc,1)),:) = Frc;
        fsize = fsize+size(Frc,1);
      end
  
      function appendV(Vrc)
        if vsize+size(Vrc,1) > size(V,1)
          Vtemp = V;
          V = [V;zeros(size(V,1),3)];
        end
        V(vsize+(1:size(Vrc,1)),:) = Vrc;
        vsize = vsize+size(Vrc,1);
      end
  
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
        if tops
          Vrc = [ ...
            X(r,1)-0.5 Y(r,1)-0.5 Z(r,1); ...
            X(r,1)-0.5 Y(r,1)+0.5 Z(r,1); ...
            X(r,1)+0.5 Y(r,1)+0.5 Z(r,1); ...
            X(r,1)+0.5 Y(r,1)-0.5 Z(r,1); ...
            ];
          Frc = [1 2 4;4 2 3];
          appendF(vsize+Frc);
          appendV(Vrc);
        end
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
      
          % Despite repeating vertices this is actually faster than adding logic
          % to reuse previous vertices 
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
      
      
          if tops
            Frc = [Frc; ...
              [1 [2 1]-2-numel(Zabove); ...
              2 1 [1]-2-numel(Zabove); ...
              ]+size(Vrc,1)];
            Vrc = [Vrc; ...
              X(r,c+1)+0.5 Y(r,c+1)-0.5 Z(r,c+1); ...
              X(r,c+1)+0.5 Y(r,c+1)+0.5 Z(r,c+1)];
          end
      
          appendF(vsize+Frc);
          appendV(Vrc);
        end
      end
      F = F(1:fsize,:);
      V = V(1:vsize,:);
    end
  
    w = size(im,2);
    h = size(im,1);
    [X,Y] = meshgrid(1:w,1:h);
    Z = im;
    [VLR,FLR] = strips(X,Y,Z,true);
    [VTB,FTB] = strips(Y',X',Z',false);
    VTB = VTB(:,[2 1 3]);
    V = [VLR;VTB];
    F = [fliplr(FLR);size(VLR,1)+FTB];
  
    % clean up
    [V,~,IM] = remove_duplicate_vertices(V,0);
    F = IM(F);
    [~,I] = unique(sort(F,2),'rows');
    F = F(I,:);
    F = F(doublearea(V,F)>0,:);
    [V,IM] = remove_unreferenced(V,F);
    F = IM(F);
  end
  function [V,F] = box_height_field_fast(im)
    % Fast version that produces boundaries between pillars

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

    [V,~,IM] = remove_duplicate_vertices(V,0);
    F = IM(F);
    [~,I] = unique(sort(F,2),'rows');
    F = F(I,:);
    F = F(doublearea(V,F)>0,:);
    [V,IM] = remove_unreferenced(V,F);
    F = IM(F);
  end


  closed = true;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Closed'},{'closed'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  if closed
    [V,F] = box_height_field_closed(im);
  else
    [V,F] = box_height_field_fast(im);
  end


end
