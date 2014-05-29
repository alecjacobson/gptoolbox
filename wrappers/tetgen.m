function [V,T,F] = tetgen(SV,SF,varargin)
  % TETGEN Call tetgen to construct a tetrahedral volume mesh with in a given
  % triangle mesh with optional internal contrained vertices.
  %
  % [V,T,F] = tetgen(SV,SF);
  % [V,T,F] = tetgen(SV,SF,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   SV  list of surface vertex positions of exterior mesh, # vertices by 3
  %   SF  list of surface face indices of exterior triangle mesh, # faces by 3
  %   Optional:
  %     'Flags'  followed by tetgen flags {'-q100'}
  % Outputs:
  %   V  list of tetrahedra vertices
  %   T  list of tetrahedra indices
  %   F  list of faces of 3D volume mesh
  %


  if ~isempty(varargin)
    assert(ischar(varargin{1}), ...
      'First optional arg not char. Using obsolete interface?');
  end

  % default values
  flags = '-q100';
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
  %off_filename = [prefix '.off'];
  %writeOFF(off_filename,SV,SF);
  
  % Try to mesh with all faces included directly

  prefix = tempname;
  poly_filename = [prefix '.poly'];
  writePOLY_tetgen(poly_filename,SV,SF,[],'BoundaryMarkers',ones(size(SF,1),1));

  %% if there are internal constraint vertices then print them to a .node file
  %if(internal_constraints)
  %  inode_filename = [prefix '.a.node'];
  %  writeNODE(inode_filename,IV);
  %end

  % graded: -q100, very-fine:-q1
  mesh_flags = '-Cpg ';
  %if(internal_constraints)
  %  flags = [flags ' -i'];
  %end
  %if(~exist('allow_resampling','var') || ~allow_resampling)
  %  flags = [flags ' -Y' '-V'];
  %end
  % call tetgen
  command = [path_to_tetgen ' ' mesh_flags ' ' flags ' ' poly_filename];
  %fprintf(command);
  [status, result] = system(command);
  if status~=0
    error(result)
  end
  % tetgen always writes output to file:
  %   xxxx.1.ele  tetrahedra
  %   xxxx.1.node tetrahedra vertices
  %   xxxx.1.face  surface faces
  ele_filename = [prefix '.1.ele'];
  face_filename = [prefix '.1.face'];
  node_filename = [prefix '.1.node'];

  F = readFACE(face_filename);
  % reverse faces because tetgen uses backwards order
  F = fliplr(F);
  % I guess this is 1-indexed because we're using a .off file rather than a
  % .poly file
  T = readELE(ele_filename);
  V = readNODE(node_filename);
  if min(T(:)) == 0 && max(T(:))<size(V,1)
    error
    % make 1-indexed
    T = T + 1;
  else if min(T(:)) >= 1
    %warning('min(T) >= 1, leaving indices as is');
  end

  delete(poly_filename);
  %if(internal_constraints)
  %  delete(inode_filename);
  %end
  delete(ele_filename);
  delete(face_filename);
  delete(node_filename);

end
