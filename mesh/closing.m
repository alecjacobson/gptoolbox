function [CV,CF,OV,OF] = closing(V,F,sigma,varargin)
  % CLOSING Create a mesh of the morphological closing of a given mesh (V,F).
  % Conducted by dilating by simga+epsilon and then eroding by sigma
  %
  % [CV,CF] = closing(V,F,sigma)
  % [CV,CF,OV,OF] = closing(V,F,sigma,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   sigma  radius of structuring element
  %   Optional:
  %     'Epsilon'  followed by epsilon value {small number % of sigma}
  %     'SignedDistanceType'  see signed_distance_isosurface.m
  %       {'pseudonormal'}
  %     'ContouringMethod'  see signed_distance_isosurface.m
  %       {'marching_cubes'}
  %     'GridSize'  see signed_distance_isosurface.m {100}
  % Outputs:
  %   CV  #CV by 3 list of closing vertex positions
  %   CF  #CF by 3 list of triangle indices into CV
  %   OV  #OV by 3 list of intermediary dilation vertex positions
  %   OF  #OF by 3 list of triangle indices into OV
  %

  

  epsilon = [];
  st = 'pseudonormal';
  ct = 'marching_cubes';
  gs = 100;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Epsilon','SignedDistanceType','ContouringMethod','GridSize'}, ...
    {'epsilon','st','ct','gs'});
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
  assert(strcmp(ct,'marching_cubes'),'alt contouring not really supported yet');

  if isempty(epsilon)
    h = max(max(V)-min(V))/gs;
    epsilon = sign(sigma)*max(0.01*sigma,0.1*h);
  end

  dim = size(V,2);
  assert(dim == size(F,2));

  switch dim
  case 2
    E = F;
    BB = bounding_box(V,'Epsilon',abs(sigma));
    [GV,side] = voxel_grid(BB,gs,'Pad',1);
    [X,Y] = deal(reshape(GV(:,1),side([2 1])),reshape(GV(:,2),side([2 1])));
    D = reshape(signed_distance([X(:) Y(:)],V,E),side([2 1]));
    [OV,OE] = contour_edges(X,Y,D,sigma);
    if isempty(OE)
      CV = [];
      CF = [];
    else
      % No real need to use the same grid. In fact the closing should stay in
      % the convex hull so we could use a tighter one.
      D = reshape(signed_distance([X(:) Y(:)],OV,OE),side([2 1]));
      [CV,CE] = contour_edges(X,Y,D,-sigma);
    end

    %surf(X,Y,0*X,'CData',D-sigma,fphong);
    %hold on;
    %tsurf(OE,OV);
    %tsurf(CE,CV);
    %tsurf(E,V);
    %hold off;
    %axis equal;
    %colormap(circshift(lipmanya(20),10,1))
    %caxis([-1 1]*max(abs(caxis)));
    

    CF = CE;
  case 3
    [OV,OF] = signed_distance_isosurface( ...
      V,F, ...
      'Level',(sigma+epsilon), ...
      'SignedDistanceType',st,'ContouringMethod','marching_cubes','GridSize',gs);
    if isempty(OF)
      CV = [];
      CF = [];
    else
      [CV,CF] = signed_distance_isosurface( ...
        OV,OF, ...
        'Level',-sigma, ...
        'SignedDistanceType',st,'ContouringMethod','marching_cubes','GridSize',gs);
    end
  end


end
