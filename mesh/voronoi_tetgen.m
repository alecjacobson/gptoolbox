function [V,E,F,C] = voronoi_tetgen(P,varargin)
  % VORONOI_TETGEN Compute a voronoi diagram using tetgen
  %
  % 
  % [V,E,F,C] = voronoi_tetgen(P,varargin)
  % 
  % Inputs:
  %   P  #P by 3 list of points
  % Outputs:
  %   V  #V by 3 list of orthocenters (?)
  %   E  #E by 2 list of edge indices into V
  %   F  #F by 3 list of face indices into ?
  %   C  #C list of cells (not supported)
  %
  % I don't trust the tetgen output...
  % 

  % default values
  flags = '';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Flags'}, {'flags'});
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

  % get a temporary file name prefix
  prefix = tempname;
  prefix = tempname;
  node_filename = [prefix '.node'];
  writeNODE(node_filename,P);

  % call tetgen
  voronoi_flags = '-vNEF';
  command = [path_to_tetgen ' ' voronoi_flags ' ' flags ' ' node_filename];
  [status, result] = system(command);
  if status~=0
    error(result)
  end

  edge_filename = [prefix '.1.v.edge'];
  face_filename = [prefix '.1.v.face'];
  node_filename = [prefix '.1.v.node'];
  cell_filename = [prefix '.1.v.cell'];

  F = readFACE(face_filename,'ForceNoBoundary',true);
  E = readEDGE(edge_filename);
  V = readNODE(node_filename);
  C = [];
  if nargout>3
    error('Cell output not handled yet.');
  end

end
