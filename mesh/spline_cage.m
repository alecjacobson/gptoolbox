function [CV,CF] = spline_cage(P,C,varargin)
  % SPLINE_CAGE Build a triangle mesh "cage" surrounding an input spline curve
  % trying to keep control points inside the triangles.
  %
  % [CV,CF] = spline_cage(P,C)
  %
  % Inputs:
  %   P  #P by 2 list of control point locations
  %   C  #C by 4 list of cubic Bezier curves as indices into rows of P
  %   Optional:
  %     'GridSize' followed by number of cells along the x-axis in the grid used
  %     for unsigned distance calculation and bwmesh {1000}
  %     'Pad' followed by number of cells to pad grid by {4}
  %     'PolyTol' followed by tolerance parameter for `spline_to_poly` used for
  %     approximating spline distances {0.1}
  %     'Tol'  followed by tolerance parameter for `bwmesh` {1.0}
  %     'SmoothingIters'  followed by smoothing iters parameter for `bwmesh`
  %       {10}
  % Outputs:
  %   CV  #CV by 2 list of cage triangulation vertex positions
  %   CF  #CF by 3 list of cage triangulation triangle indices into rows of CV
  %

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'GridSize','Pad','PolyTol','Tol','SmoothingIters'}, ...
    {'gs','pad','poly_tol','bwmesh_tol','bwmesh_smoothing_iters'});
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

  gs = 1000;
  pad = 4;
  poly_tol = 0.1;
  bwmesh_tol = 1.0;
  bwmesh_smoothing_iters = 10;

  % Just use distance to a well sampled poly-line
  [V,E] = spline_to_poly(P,C,poly_tol);
  [BC,side,r] = voxel_grid(P,gs,'Pad',pad);
  X = reshape(BC(:,1),side([2 1]));
  Y = reshape(BC(:,2),side([2 1]));
  sqrD = point_mesh_squared_distance(BC,V,E);
  D = sqrt(reshape(sqrD,size(X)));
  % Find maximum distance sampled at control points to try to be sure they're
  % inside the extracted mesh
  Dmax = max(interp2(X,Y,reshape(D,size(X)),P(:,1),P(:,2)));
  % convert into fuzzy logistic looking function
  h = 2/sqrt(2)*r(1);
  Dmax = Dmax+h;
  Z = 1-min(max(D-Dmax,0),h)/h;
  % Ignore holes
  Z = imfill(Z);
  [CV,CF] = bwmesh(Z,'Tol',bwmesh_tol,'SmoothingIters',bwmesh_smoothing_iters);
  % bwmesh flips things around and doesn't know the scale
  CV(:,2) = size(Z,1)-CV(:,2);
  CV = (CV-0.0)/size(Z,2)*(max(X(:))-min(X(:)))+[min(X(:)) min(Y(:))];
end
