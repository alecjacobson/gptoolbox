function [DV,DT,DF,b,bc] = encage(V,CV,CF,varargin)
  % ENCAGE  Create a mesh and compute boundary conditions for computing weights
  % on for points V embedded in a polyhedron cage (CV,CF).
  %
  % [DV,DT,DF,b,bc] = encage(V,CV,CF)
  % [DV,DT,DF,b,bc] = encage(V,CV,CF,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   CV  #CV by dim list of cage vertex positions
  %   CF  #CF by dim list of cage face indices into CV
  %   Optional:
  %     'TetgenFlags' followed by a string of tetgen flags {'-q100'}
  %     'Quiet' followed by whether to tell tetgen to be quiet {false}
  % Outputs:
  %   DV  #DV by dim list of tet mesh vertex positions
  %   DT  #DT by dim+1 list of tet mesh tet indices into DV
  %   DF  #DF by dim list of tet mesh boundary facet indices into DV
  %   b  #b list of indices of boundary vertices in V
  %   bc  #b by #CV list of boundary conditions
  % 

  tetgen_flags = '-q100';
  quiet = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'TetgenFlags','Quiet'}, ...
    {'tetgen_flags','quiet'});
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

  if quiet
    fid = fopen('/dev/null','w');
  else
    fid = 1;
  end

  prefix = tempname;
  poly_filename = [prefix '.poly'];
  writePOLY_tetgen(poly_filename,[CV;V],CF,[],'BoundaryMarkers',(1:size(CF,1))');

  delaunay_flags = '-p';
  command = ...
    [path_to_tetgen ' ' delaunay_flags ' ' tetgen_flags ' ' poly_filename];
  fprintf(fid,'%s\n',command);
  if quiet
    [status, result] = system(command);
  else
    [status, result] = system(command,'-echo');
  end
  fprintf(fid,'%s',result);
  if status ~= 0
    fprintf(fid,'%s',result);
  else
    if strfind(result,'Jettisoning redundant points.')
      fprintf(fid,'%s',result);
      warning('Tetgen removed duplicate points');
    end
  end
  switch(status)
  case 0
  case 3
    error('Self intersections. Try selfintersect() or remove self-intersections');
  case 4
  % Small feature
    error('Try again with reduced T?');
  case 5
  % close input facets
    error('Try again with -Y');
  case 6
  % close input facets
    error('Input error. Check yo` input');
  case 134
    error('Non-planar facets?');
  otherwise
    error('Tetgen returned status %d != 0',status);
  end
  delete(poly_filename);

  ele_filename = [prefix '.1.ele'];
  face_filename = [prefix '.1.face'];
  node_filename = [prefix '.1.node'];

  DV = readNODE(node_filename);
  DT = readELE(ele_filename);
  [DF,BF] = readFACE(face_filename);

  delete(ele_filename);
  delete(node_filename);
  delete(face_filename);

  assert(size(DF,1) == size(BF,1));
  % B2CF(i,:) = [v f] --> vertex v in DV is on face f in CF
  B2CF = [DF(:) repmat(BF,size(DF,2),1)];
  % Just keep one matching for each vertex
  [~,UI] = unique(B2CF(:,1));
  B2CF = B2CF(UI,:);
  % boundary conditions
  b = B2CF(:,1);
  if size(CF,2) == 3
    % Corners of face for each vertex
    K = CF(B2CF(:,2),:);
    B = barycentric_coordinates( ...
      DV(B2CF(:,1),:), CV(K(:,1),:), CV(K(:,2),:), CV(K(:,3),:));
    bc = sparse(repmat(1:numel(b),1,3)',K(:),B(:),numel(b),size(CV,1));
  else
    if nargout>4
      error('bc not supported for non-triangle mesh cages');
    end
  end
  % Faces are flipped
  DF = fliplr(DF);

end
