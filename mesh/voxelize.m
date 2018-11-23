function [W,BC,DV,Q,r] = voxelize(V,F,side,varargin)
  % VOXELIZE Given a mesh compute a voxelization of that mesh on a regular grid
  % fit to the bounding box.
  %
  % W = voxelize(V,F,side)
  % [W,BC,DV,Q] = voxelize(V,F,side,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into V, expects counterclockwise
  %     orientation to produce outward pointing normal
  %   side  either the number of cells in the x direction or a triple
  %     containing the number of cells in each dimension
  %   Optional:
  %     'Pad'  followed by number of padding to add to each side. 1 would mean
  %       add one extra cell on the top, bottom, left, right, front, and back
  %       of bounding box lattice {0}.
  %     'Boundary'  followed by wether to consider any cell intersecting the
  %       mesh as inside {true}.
  %     'Interior'  followed by whether to consider any cell full inside the
  %        mesh as inside {true}
  %     'Closed'  followed by whether to assume mesh is closed when determining
  %       interior (avoid winding number computation, use flood filling)
  %     'Centers'  followed by prod(side) centers (e.g., already produced by
  %       call to previous call to voxel_grid.m)
  % Outputs:
  %   W  side(1) by side(2) by side(3) matrix with W(i,j,k) ~= if location
  %     cell centered at BC(i,j,k) overlaps with the volume of (V,F)
  %   BC  prod(side) by dim matrix of cell barycenters
  %   DV  side(1)+1 by side(2)+1 by side(3)+1 matrix of cell corner locations
  %   Q  #Q by 4 list of quads indexing DV
  %   r  dim-long list of voxel widths in each direction
  %
  % Known issues: the ouput surface mesh will contain non-manifold edges and
  % vertices at voxels that meet at such edges or vertices. 
  %
  % Example:
  %   [W,BC,DV,Q] = voxelize(V,F,50);
  %   % Get triangles from quads
  %   DF = [Q(:,[1 2 3]);Q(:,[1 3 4])];
  %   % Remove unreferenced corners
  %   [SV,IM] = remove_unreferenced(DV,DF);
  %   % re-index
  %   SF = IM(DF);
  %   SQ = IM(Q);
  %
  %   [W,BC] = voxelize(V,F,50,'Pad',1);
  %   % Use matlab's isosurface to extract surface as triangle mesh
  %   BC = reshape(BC,[size(W) 3]);
  %   surf = isosurface(BC(:,:,:,1),BC(:,:,:,2),BC(:,:,:,3),W,0.5);
  %   DV = bsxfun(@plus,bsxfun(@times,bsxfun(@rdivide,surf.vertices-1, ...
  %     [size(W,2) size(W,1) size(W,3)]-1),max(BC)-min(BC)),min(BC));
  %   DF = surf.faces;
  % 
  %
  % See also: bwboundaries, isosurface

  with_boundary = true;
  with_interior = true;
  pad_count = 0;
  closed = [];
  BC = [];
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Boundary','Interior','Pad','Centers','Closed'}, ...
    {'with_boundary','with_interior','pad_count','BC','closed'});
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
  if isempty(closed)
    closed = isempty(boundary_faces(F));
  end


  assert(any(size(side) == 1),'side should be a vector');
  % Make sure we have a row vector
  side = side(:)';

  dim = size(V,2);

  if isempty(BC)
    [BC,side,r] = voxel_grid(V,side,'Pad',pad_count);
  else
    assert(prod(side) == size(BC,1),'side dims must match BC');
    r = (max(BC)-min(BC))./(side-1);
  end
  NV = min(BC);
  XV = max(BC);

  switch dim
  case 2
    W  = zeros([side(2) side(1)]);
    [DV,~,QT,QL] = voxel_surface(W,'Centers',BC);

    if with_interior && ~closed
      WDV = reshape(winding_number(V,F,DV),[side(2) side(1)]+1);
      WDV = abs(WDV) >= 0.5;
      W = ( ...
        WDV(1:end-1,1:end-1) | ...
        WDV(1:end-1,  2:end) | ...
        WDV(  2:end,1:end-1) | ...
        WDV(  2:end,  2:end));
    end
    Q = [QT;QL];
    if with_boundary
      for ei = 1:size(F,1)

        sqrD = segment_segment_squared_distance( ...
          V(F(ei,1),:), V(F(ei,2),:), ...
          DV(Q(:,1),:), DV(Q(:,2),:));
        IF = mod(find(sqrD<eps)-1,size(Q,1))+1;
        iQT = IF(IF<=size(QT,1));
        iQL = IF(IF>size(QT,1))-size(QT,1);

        [II,JJ] = ind2sub([side(2) side(1)]+1,max(QT(iQT,:),[],2));
        W(sub2ind(side([2 1]),II-1,JJ-1)) = 1;
        W(sub2ind(side([2 1]),II  ,JJ-1)) = 1;

        [II,JJ] = ind2sub([side(2) side(1)]+1,max(QL(iQL,:),[],2));
        W(sub2ind(side([2 1]),II-1,JJ-1)) = 1;
        W(sub2ind(side([2 1]),II-1,JJ  )) = 1;

        %clf;
        %hold on;
        %surf(reshape(BC(:,1),size(W)),reshape(BC(:,2),size(W)),reshape(BC(:,1)*0,size(W)),'CData',W*1,fphong,'EdgeColor','none');
        %plot_edges(V,F(ei,:),'y','LineWidth',8);
        %plot_edges(V,F,'g','LineWidth',4);
        %text(BC(:,1),BC(:,2),num2str((1:size(BC,1))'),'BackgroundCOlor',[0.8 0.8 0.8]);
        %text(DV(:,1),DV(:,2),num2str((1:size(DV,1))'),'BackgroundCOlor',[0.8 0.1 0.1]);
        %plot_edges(DV,QT,'--r','LineWidth',1);
        %plot_edges(DV,QL,'--b','LineWidth',1);
        %plot_edges(DV,QT(iQT,:),'r','LineWidth',3);
        %plot_edges(DV,QL(iQL,:),'b','LineWidth',3);
        %hold off;
        %axis equal;
        %colormap([zeros(20,3);1 1 1]);
        %drawnow;

      end
    end
    if with_interior && closed
      W = imfill(W,'holes');
    end
    % Pad winding number with 0s and take differences in all 6 directions
    Wp = padarray(W,[1 1],0);
    Dy = diff(Wp,1,1);
    Dx = diff(Wp,1,2);
    Q = [ ...
      QL(Dx(2:end-1,1:end,2:end-1)>0.5,:); ...
      QT(Dy(1:end,2:end-1,2:end-1)>0.5,:); ...
      fliplr(QL(Dx(2:end-1,1:end,2:end-1)<-0.5,:)); ...
      fliplr(QT(Dy(1:end,2:end-1,2:end-1)<-0.5,:)); ...
      ];
  case 3
    W = zeros([side(2) side(1) side(3)]);
    [DV,~,QT,QF,QL] = voxel_surface(W,'Centers',BC);
  
    if with_interior && ~closed
      % Winding number is a heavy handed way of determining inside/outside for a
      % simple closed polyhedron.  If the boundary has already been detected then
      % this should be floodfilling instead.
  
      %W = winding_number(V,F,BC);
      WDV = reshape(winding_number(V,F,DV),[side(2) side(1) side(3)]+1);
      WDV = abs(WDV) >= 0.5;
      W = ( ...
        WDV(1:end-1,1:end-1,1:end-1) | ...
        WDV(1:end-1,1:end-1,2:end) | ...
        WDV(1:end-1,  2:end,1:end-1) | ...
        WDV(1:end-1,  2:end,2:end) | ...
        WDV(  2:end,1:end-1,1:end-1) | ...
        WDV(  2:end,1:end-1,2:end) | ...
        WDV(  2:end,  2:end,1:end-1) | ...
        WDV(  2:end,  2:end,2:end));
    end
  
    Q = [QT;QF;QL];
      
    if with_boundary
      % TODO: Could at least ignore already-inside cells.
      DF = [Q(:,[3 2 1]);Q(:,[4 3 1])];
      IF = intersect_other(V,F,DV,DF);
      IF = mod(IF(:,2)-1,size(Q,1))+1;
      W = reshape(W,[side(2) side(1) side(3)]);
      % Force winding number to mark voxels intersecting boundary as inside
      % (This should just be replaced with rasterizing the surface)
      iQT = IF(IF<=size(QT,1));
      iQF = IF(IF>size(QT,1) & IF<=size(QT,1)+size(QF,1))-size(QT,1);
      iQL = IF(IF>size(QT,1)+size(QF,1))-size(QT,1)-size(QF,1);
      %W(:) = 0;
      %warning('winding number = 0');
  
      [II,JJ,KK] = ind2sub([side(2)+1 side(1)   side(3)  ],iQT);
      KK = KK(II<=side(2));
      JJ = JJ(II<=side(2));
      II = II(II<=side(2));
      W(sub2ind([side(2) side(1) side(3)],II,JJ,KK)) = 1;
      [II,JJ,KK] = ind2sub([side(2)+1 side(1)   side(3)  ],iQT);
      II = II-1;
      KK = KK(II<=side(2));
      JJ = JJ(II<=side(2));
      II = II(II<=side(2));
      W(sub2ind([side(2) side(1) side(3)],II,JJ,KK)) = 1;
  
      [II,JJ,KK] = ind2sub([side(2)   side(1)+1 side(3)  ],iQL);
      II = II(JJ<=side(1));
      KK = KK(JJ<=side(1));
      JJ = JJ(JJ<=side(1));
      W(sub2ind([side(2) side(1) side(3)],II,JJ,KK)) = 1;
      [II,JJ,KK] = ind2sub([side(2)   side(1)+1 side(3)  ],iQL);
      JJ = JJ-1;
      II = II(JJ<=side(1));
      KK = KK(JJ<=side(1));
      JJ = JJ(JJ<=side(1));
      W(sub2ind([side(2) side(1) side(3)],II,JJ,KK)) = 1;
  
      [II,JJ,KK] = ind2sub([side(2)   side(1)   side(3)+1],iQF);
      II = II(KK<=side(3));
      JJ = JJ(KK<=side(3));
      KK = KK(KK<=side(3));
      W(sub2ind([side(2) side(1) side(3)],II,JJ,KK)) = 1;
      [II,JJ,KK] = ind2sub([side(2)   side(1)   side(3)+1],iQF);
      KK = KK-1;
      II = II(KK<=side(3));
      JJ = JJ(KK<=side(3));
      KK = KK(KK<=side(3));
      W(sub2ind([side(2) side(1) side(3)],II,JJ,KK)) = 1;
    end
  
    if with_interior && closed
      W = imfill(W,'holes');
    end
  
    %trisurf(Q(IF,:),DV(:,1),DV(:,2),DV(:,3),'FaceColor','r');
    %hold on;
    %tsurf(F,V,'FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.2);
    %hold off;
  
    %% From cells to faces:
    %[II,JJ,KK] = ind2sub([side(2) side(1) side(3)],find(abs(W)>0.5));
    %Q = [ ...
    %  QT(sub2ind([side(2)+1 side(1)   side(3)  ],II+1,JJ,KK),:); ...
    %  QT(sub2ind([side(2)+1 side(1)   side(3)  ],II+0,JJ,KK),:); ...
    %  QF(sub2ind([side(2)   side(1)   side(3)+1],II,JJ,KK+1),:); ...
    %  QF(sub2ind([side(2)   side(1)   side(3)+1],II,JJ,KK+0),:); ...
    %  QL(sub2ind([side(2)   side(1)+1 side(3)  ],II,JJ+1,KK),:); ...
    %  QL(sub2ind([side(2)   side(1)+1 side(3)  ],II,JJ+0,KK),:); ...
    %  ];
    %trisurf(Q,DV(:,1),DV(:,2),DV(:,3),'FaceAlpha',0.5,'EdgeAlpha',0.5);
    %hold on;
    %tsurf(F,V,'FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.2);
    %hold off;
    %axis equal
    %
    %error
  
    % Pad winding number with 0s and take differences in all 6 directions
    Wp = padarray(W,[1 1 1],0);
    Dy = diff(Wp,1,1);
    Dx = diff(Wp,1,2);
    Dz = diff(Wp,1,3);
    Q = [ ...
      QF(Dz(2:end-1,2:end-1,1:end)>0.5,:); ...
      QL(Dx(2:end-1,1:end,2:end-1)>0.5,:); ...
      QT(Dy(1:end,2:end-1,2:end-1)>0.5,:); ...
      fliplr(QF(Dz(2:end-1,2:end-1,1:end)<-0.5,:)); ...
      fliplr(QL(Dx(2:end-1,1:end,2:end-1)<-0.5,:)); ...
      fliplr(QT(Dy(1:end,2:end-1,2:end-1)<-0.5,:)); ...
      ];
  end
end
