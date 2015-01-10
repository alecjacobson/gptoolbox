function [W,BC,DV,Q] = voxelize(V,F,side,varargin)
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
  % Outputs:
  %   W  side(1) by side(2) by side(3) matrix with W(i,j,k) ~= if location
  %     cell centered at BC(i,j,k) overlaps with the volume of (V,F)
  %   BC  side(1) by side(2) by side(3) matrix of cell barycenters
  %   DV  side(1)+1 by side(2)+1 by side(3)+1 matrix of cell corner locations
  %   Q  #Q by 4 list of quads indexing DV
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
  %   surf = isosurface(BC3(:,:,:,1),BC3(:,:,:,2),BC3(:,:,:,3),W,0.5);
  %   DF = surf.faces;
  %   DV = surf.vertices;
  % 
  %
  % See also: bwboundaries, isosurface

  with_boundary = true;
  with_interior = true;
  pad_count = 0;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Boundary','Interior','Pad'}, ...
    {'with_boundary','with_interior','pad_count'});
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


  NV = min(V);
  XV = max(V);

  switch numel(side)
  case 3
    side_x = side(1);
    side_y = side(2);
    side_z = side(3);
  case 1
    % Enclose bounding box in regular mesh
    side_x = side(1);
    side_y = round(side_x * (XV(2)-NV(2))/(XV(1)-NV(1)));
    side_z = round(side_x * (XV(3)-NV(3))/(XV(1)-NV(1)));
  otherwise
    error('side should be scalar or triplet');
  end
  assert(all(side==ceil(side)),'All side values should be integer');

  pad_count = 0;
  pad_count = 1;
  pad_length = pad_count*(XV(1)-NV(1))/(side_x-2*pad_count);
  NV = NV-pad_length;
  XV = XV+pad_length;

  [X,Y,Z] = meshgrid( ...
    NV(1)+linspace(0,1,side_x)*(XV(1)-NV(1)), ...
    NV(2)+linspace(0,1,side_y)*(XV(2)-NV(2)), ...
    NV(3)+linspace(0,1,side_z)*(XV(3)-NV(3)));
  % indices of cells
  I = reshape(1:numel(X),[side_y side_x side_z]);
  % barycenters of cells
  BC = [X(:) Y(:) Z(:)];
  r = [max(V)-min(V)]./([side_x side_y side_z] - 1);
  [X,Y,Z] = meshgrid( ...
    (NV(1)-0.5*r(1))+linspace(0,1,side_x+1)*(XV(1)-NV(1)+r(1)), ...
    (NV(2)-0.5*r(2))+linspace(0,1,side_y+1)*(XV(2)-NV(2)+r(2)), ...
    (NV(3)-0.5*r(3))+linspace(0,1,side_z+1)*(XV(3)-NV(3)+r(3)));
  % corners of cells
  DV = [X(:) Y(:) Z(:)];

  if with_interior
    % Winding number is a heavy handed way of determining inside/outside for a
    % simple closed polyhedron.  If the boundary has already been detected then
    % this should be floodfilling instead.

    %W = winding_number(V,F,BC);
    WDV = reshape(winding_number(V,F,DV),[side_y side_x side_z]+1);
    WDV = WDV >= 0.5;
    W = ( ...
      WDV(1:end-1,1:end-1,1:end-1) | ...
      WDV(1:end-1,1:end-1,2:end) | ...
      WDV(1:end-1,  2:end,1:end-1) | ...
      WDV(1:end-1,  2:end,2:end) | ...
      WDV(  2:end,1:end-1,1:end-1) | ...
      WDV(  2:end,1:end-1,2:end) | ...
      WDV(  2:end,  2:end,1:end-1) | ...
      WDV(  2:end,  2:end,2:end));
  else
    W = zeros(size(BC));
  end

  if with_boundary
    % This is a very stupid way to implement this. The correct thing to do is
    % to _rasterize_ each triangle. Instead this is building a list of all
    % interfaces in the voxel mesh and literally checking intersections against
    % the input mesh.
    [II,JJ,KK] = ind2sub([side_y side_x side_z]+1,reshape(1:size(DV,1),[side_y side_x side_z]+1));

    QF = [ ...
      reshape(sub2ind([side_y side_x side_z]+1,II(1:end-1,1:end-1,1:end),JJ(1:end-1,1:end-1,1:end),KK(1:end-1,1:end-1,1:end)),[],1) ...
      reshape(sub2ind([side_y side_x side_z]+1,II(  2:end,1:end-1,1:end),JJ(  2:end,1:end-1,1:end),KK(  2:end,1:end-1,1:end)),[],1) ...
      reshape(sub2ind([side_y side_x side_z]+1,II(  2:end,  2:end,1:end),JJ(  2:end,  2:end,1:end),KK(  2:end,  2:end,1:end)),[],1) ...
      reshape(sub2ind([side_y side_x side_z]+1,II(1:end-1,  2:end,1:end),JJ(1:end-1,  2:end,1:end),KK(1:end-1,  2:end,1:end)),[],1) ...
      ];

    QL = fliplr([ ...
      reshape(sub2ind([side_y side_x side_z]+1,II(1:end-1,1:end,1:end-1),JJ(1:end-1,1:end,1:end-1),KK(1:end-1,1:end,1:end-1)),[],1) ...
      reshape(sub2ind([side_y side_x side_z]+1,II(  2:end,1:end,1:end-1),JJ(  2:end,1:end,1:end-1),KK(  2:end,1:end,1:end-1)),[],1) ...
      reshape(sub2ind([side_y side_x side_z]+1,II(  2:end,1:end,  2:end),JJ(  2:end,1:end,  2:end),KK(  2:end,1:end,  2:end)),[],1) ...
      reshape(sub2ind([side_y side_x side_z]+1,II(1:end-1,1:end,  2:end),JJ(1:end-1,1:end,  2:end),KK(1:end-1,1:end,  2:end)),[],1) ...
      ]);

    QT = [ ...
      reshape(sub2ind([side_y side_x side_z]+1,II(1:end,1:end-1,1:end-1),JJ(1:end,1:end-1,1:end-1),KK(1:end,1:end-1,1:end-1)),[],1) ...
      reshape(sub2ind([side_y side_x side_z]+1,II(1:end,  2:end,1:end-1),JJ(1:end,  2:end,1:end-1),KK(1:end,  2:end,1:end-1)),[],1) ...
      reshape(sub2ind([side_y side_x side_z]+1,II(1:end,  2:end,  2:end),JJ(1:end,  2:end,  2:end),KK(1:end,  2:end,  2:end)),[],1) ...
      reshape(sub2ind([side_y side_x side_z]+1,II(1:end,1:end-1,  2:end),JJ(1:end,1:end-1,  2:end),KK(1:end,1:end-1,  2:end)),[],1) ...
      ];
    Q = [QT;QF;QL];
    % TODO: Could at least ignore already-inside cells.
    DF = [Q(:,[3 2 1]);Q(:,[4 3 1])];
    IF = intersect_other(V,F,DV,DF);
    IF = mod(IF(:,2)-1,size(Q,1))+1;
    W = reshape(W,[side_y side_x side_z]);
    % Force winding number to mark voxels intersecting boundary as inside
    % (This should just be replaced with rasterizing the surface)
    iQT = IF(IF<=size(QT,1));
    iQF = IF(IF>size(QT,1) & IF<=size(QT,1)+size(QF,1))-size(QT,1);
    iQL = IF(IF>size(QT,1)+size(QF,1))-size(QT,1)-size(QF,1);
    %W(:) = 0;
    %warning('winding number = 0');

    [II,JJ,KK] = ind2sub([side_y+1 side_x   side_z  ],iQT);
    KK = KK(II<=side_y);
    JJ = JJ(II<=side_y);
    II = II(II<=side_y);
    W(sub2ind([side_y side_x side_z],II,JJ,KK)) = 1;
    [II,JJ,KK] = ind2sub([side_y+1 side_x   side_z  ],iQT);
    II = II-1;
    KK = KK(II<=side_y);
    JJ = JJ(II<=side_y);
    II = II(II<=side_y);
    W(sub2ind([side_y side_x side_z],II,JJ,KK)) = 1;

    [II,JJ,KK] = ind2sub([side_y   side_x+1 side_z  ],iQL);
    II = II(JJ<=side_x);
    KK = KK(JJ<=side_x);
    JJ = JJ(JJ<=side_x);
    W(sub2ind([side_y side_x side_z],II,JJ,KK)) = 1;
    [II,JJ,KK] = ind2sub([side_y   side_x+1 side_z  ],iQL);
    JJ = JJ-1;
    II = II(JJ<=side_x);
    KK = KK(JJ<=side_x);
    JJ = JJ(JJ<=side_x);
    W(sub2ind([side_y side_x side_z],II,JJ,KK)) = 1;

    [II,JJ,KK] = ind2sub([side_y   side_x   side_z+1],iQF);
    II = II(KK<=side_z);
    JJ = JJ(KK<=side_z);
    KK = KK(KK<=side_z);
    W(sub2ind([side_y side_x side_z],II,JJ,KK)) = 1;
    [II,JJ,KK] = ind2sub([side_y   side_x   side_z+1],iQF);
    KK = KK-1;
    II = II(KK<=side_z);
    JJ = JJ(KK<=side_z);
    KK = KK(KK<=side_z);
    W(sub2ind([side_y side_x side_z],II,JJ,KK)) = 1;
  end

  %trisurf(Q(IF,:),DV(:,1),DV(:,2),DV(:,3),'FaceColor','r');
  %hold on;
  %tsurf(F,V,'FaceColor','g','FaceAlpha',0.2,'EdgeAlpha',0.2);
  %hold off;

  %% From cells to faces:
  %[II,JJ,KK] = ind2sub([side_y side_x side_z],find(abs(W)>0.5));
  %Q = [ ...
  %  QT(sub2ind([side_y+1 side_x   side_z  ],II+1,JJ,KK),:); ...
  %  QT(sub2ind([side_y+1 side_x   side_z  ],II+0,JJ,KK),:); ...
  %  QF(sub2ind([side_y   side_x   side_z+1],II,JJ,KK+1),:); ...
  %  QF(sub2ind([side_y   side_x   side_z+1],II,JJ,KK+0),:); ...
  %  QL(sub2ind([side_y   side_x+1 side_z  ],II,JJ+1,KK),:); ...
  %  QL(sub2ind([side_y   side_x+1 side_z  ],II,JJ+0,KK),:); ...
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
