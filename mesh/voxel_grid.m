function [BC,side,r] = voxel_grid(V,side,varargin)
  % VOXEL_GRID Prepare a voxel grid around a set of points V
  % 
  % [BC,side,r] = voxel_grid(V,side,varargin)
  % 
  % Inputs:
  %   V  #V by 3 list of input point positions
  %   side  either:
  %     scalar specifying how many steps along the x-coordinate of V's bounding
  %     box
  %          or
  %     3 list specifying how many steps in x, y, and z coordinates of V's
  %       bounding box
  %   Optional:
  %     'Pad' followed by the number of "extra cells" to add on all six sides
  %       of the grid (pad_count).
  % Outputs:
  %   BC  prod(side+2*pad_count) by 3 list of cell centers
  %   side  number of cells on each side: side+2*pad_count
  %   r  size of step in each direciton
  %
  % See also: voxelize


  pad_count = 0;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Pad'}, ...
    {'pad_count'});
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

  dim = size(V,2);

  assert(all(side>(pad_count*2+1)),'side should be > 2*pad_count+1');
  side = side-pad_count*2;
  switch numel(side)
  case 3
    side(1) = side(1);
    side(2) = side(2);
    side(3) = side(3);
    NV = min(V);
    XV = max(V);
    r = [XV-NV]./([side(1) side(2) side(3)] - 1);
  case 2
    assert(size(V,2) == 2);
    side(1) = side(1);
    side(2) = side(2);
    NV = min(V);
    XV = max(V);
    r = [XV-NV]./([side(1) side(2)] - 1);
  case 1
    % Enclose bounding box in regular mesh
    side(1) = side(1);
    NV = min(V);
    XV = max(V);
    for d = 2:dim
      side(d) = ceil(side(1) * (XV(d)-NV(d))/(XV(1)-NV(1)));
    end
    r = max((XV-NV)./(side-1));
    % recenter
    old_cen = 0.5*(XV + NV);
    XV = NV + r*(side-1);
    cen = 0.5*(XV + NV);
    XV = XV + old_cen - cen;
    NV = NV + old_cen - cen;
  otherwise
    error('side should be scalar or triplet');
  end
  assert(all(side==ceil(side)),'All side values should be integer');

  side = side+pad_count*2;
  old_cen = 0.5*(XV + NV);
  XV = NV+r.*(side-1);
  cen = 0.5*(XV + NV);
  XV = XV + old_cen - cen;
  NV = NV + old_cen - cen;

  r = (XV-NV)./(side-1);

  switch dim
  case 3
    [X,Y,Z] = meshgrid( ...
      NV(1)+linspace(0,1,side(1))*(XV(1)-NV(1)), ...
      NV(2)+linspace(0,1,side(2))*(XV(2)-NV(2)), ...
      NV(3)+linspace(0,1,side(3))*(XV(3)-NV(3)));
    % barycenters of cells
    BC = [X(:) Y(:) Z(:)];
  case 2
    [X,Y] = meshgrid( ...
      NV(1)+linspace(0,1,side(1))*(XV(1)-NV(1)), ...
      NV(2)+linspace(0,1,side(2))*(XV(2)-NV(2)));
    % barycenters of cells
    BC = [X(:) Y(:)];
  end

end
