function [V,Q,QT,QF,QL] = voxel_surface(W,varargin)
  % VOXEL_SURFACE Compute the surface quad mesh of a 3d logical image.
  %
  % [V,Q,QT,QF,QL] = voxel_surface(W,varargin)
  %
  % Inputs:
  %   W  n by m by k binary matix
  %   Optional:
  %     'Centers' followed by n*m*k by 3 list of cell centers
  % Outputs:
  %   V  (n+1)*(m+1)*(k+1) by 3 list of cell corner positions
  %   Q  #Q by 4 list of quads indexing V
  %   QT  n*m*k by 4 list of "top" quad faces indexing V
  %   QF  n*m*k by 4 list of "front" quad faces indexing V
  %   QL  n*m*k by 4 list of "left" quad faces indexing V
  %
  % Examples:
  %   [V,Q] = voxel_surface(W);
  %   [V,IM] = remove_unreferenced(V,Q);
  %   Q = IM(Q);
  %   trisurf(Q,V(:,1),V(:,2),V(:,3));
  %   axis equal;
  %
  % See also: voxelize, voxel_grid

  BC = [];
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Centers'}, ...
    {'BC'});
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

  dim = ndims(W);

  switch dim
  case 2
    side = [size(W,2) size(W,1)];
    if isempty(BC)
      [X,Y] = meshgrid( ...
        linspace(0,1,side(1)), ...
        linspace(0,1,side(2)));
      % barycenters of cells
      BC = [X(:) Y(:)];
    end
    NV = min(BC);
    XV = max(BC);
    r = (XV-NV)./(side-1);
    [X,Y] = meshgrid( ...
      (NV(1)-0.5*r(1))+linspace(0,1,side(1)+1)*(XV(1)-NV(1)+r(1)), ...
      (NV(2)-0.5*r(2))+linspace(0,1,side(2)+1)*(XV(2)-NV(2)+r(2)));

    [Q,V] = surf2patch(X,Y,0*X);
    V = V(:,1:2);
    QT = [Q(:,[3 4]);Q(:,[1 2]);];
    QF = [Q(:,[2 3]); Q(:,[4 1])];
    if ~any(W(:))
      Q = [];
    else
      assert(false,'non trivial W not supported');
    end

  case 3
    side = [size(W,2) size(W,1) size(W,3)];
    if isempty(BC)
      [X,Y,Z] = meshgrid( ...
        linspace(0,1,side(1)), ...
        linspace(0,1,side(2)), ...
        linspace(0,1,side(3)));
      % barycenters of cells
      BC = [X(:) Y(:) Z(:)];
    end

    NV = min(BC);
    XV = max(BC);
    r = (XV-NV)./(side-1);
    [X,Y,Z] = meshgrid( ...
      (NV(1)-0.5*r(1))+linspace(0,1,side(1)+1)*(XV(1)-NV(1)+r(1)), ...
      (NV(2)-0.5*r(2))+linspace(0,1,side(2)+1)*(XV(2)-NV(2)+r(2)), ...
      (NV(3)-0.5*r(3))+linspace(0,1,side(3)+1)*(XV(3)-NV(3)+r(3)));
    % corners of cells
    V = [X(:) Y(:) Z(:)];

    [II,JJ,KK] = ind2sub([side(2) side(1) side(3)]+1,reshape(1:size(V,1),[side(2) side(1) side(3)]+1));

    QF = [ ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(1:end-1,1:end-1,1:end),JJ(1:end-1,1:end-1,1:end),KK(1:end-1,1:end-1,1:end)),[],1) ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(  2:end,1:end-1,1:end),JJ(  2:end,1:end-1,1:end),KK(  2:end,1:end-1,1:end)),[],1) ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(  2:end,  2:end,1:end),JJ(  2:end,  2:end,1:end),KK(  2:end,  2:end,1:end)),[],1) ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(1:end-1,  2:end,1:end),JJ(1:end-1,  2:end,1:end),KK(1:end-1,  2:end,1:end)),[],1) ...
        ];

    QL = fliplr([ ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(1:end-1,1:end,1:end-1),JJ(1:end-1,1:end,1:end-1),KK(1:end-1,1:end,1:end-1)),[],1) ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(  2:end,1:end,1:end-1),JJ(  2:end,1:end,1:end-1),KK(  2:end,1:end,1:end-1)),[],1) ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(  2:end,1:end,  2:end),JJ(  2:end,1:end,  2:end),KK(  2:end,1:end,  2:end)),[],1) ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(1:end-1,1:end,  2:end),JJ(1:end-1,1:end,  2:end),KK(1:end-1,1:end,  2:end)),[],1) ...
        ]);

    QT = [ ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(1:end,1:end-1,1:end-1),JJ(1:end,1:end-1,1:end-1),KK(1:end,1:end-1,1:end-1)),[],1) ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(1:end,  2:end,1:end-1),JJ(1:end,  2:end,1:end-1),KK(1:end,  2:end,1:end-1)),[],1) ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(1:end,  2:end,  2:end),JJ(1:end,  2:end,  2:end),KK(1:end,  2:end,  2:end)),[],1) ...
        reshape(sub2ind([side(2) side(1) side(3)]+1,II(1:end,1:end-1,  2:end),JJ(1:end,1:end-1,  2:end),KK(1:end,1:end-1,  2:end)),[],1) ...
        ];

    if ~any(W(:))
      Q = [];
    else
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
end
