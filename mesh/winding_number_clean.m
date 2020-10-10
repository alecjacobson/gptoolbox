function [CV,CF,TV,TT,w,SV,SF] = winding_number_clean(V,F,varargin)
  % WINDING_NUMBER_CLEAN Clean a mesh using the generalized winding number
  % pipeline with some default settings: "Shoot for the moon!!!!"
  %
  % [CV,CT] = winding_number_clean(V,F)
  % [CV,CT] = winding_number_clean(V,F,'ParameterName',ParameterValue,..)
  %
  % Inputs:
  %   V  #V by 3 input mesh positions
  %   F  #F by 3 list of triangle indices into V
  %   Optional:
  %     'SurfaceOnly' followed by whether to return CF  #CF by 3 list of
  %       surface triangle indices instead of CT, and remove rows unrefernced
  %       by CF from CV. Equivalent to running:
  %           [CV,CT] = winding_number_clean(V,F);
  %           CF = boundary_faces(CT);           
  %           [CV,I] = remove_unreferenced(TV,CF)
  %           CF = I(CF);
  %       {false}.
  % Outputs:
  %   CV  #CV by 3 output mesh positions
  %   CT  #CT by 3 list of tet indices into CV
  %   TV  #TV by 3 list of tetrahedral mesh vertex positions
  %   TT  #TT by 3 list of tetrahedral mesh tet indices
  %   w  #TT list of winding numbers
  %   SV  #SV by 3 list of self-intersected "clean", but non-solid mesh indices
  %   SF  #SF by 3 list of triangle indices into SV
  % 

  surface_only = false;
  % parse optional input parameters
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'SurfaceOnly'},{'surface_only'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin))
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace 
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  % shoot for the moon with winding number on default settings
  [SV,SF] = clean_mesh(V,F,'Quiet',true);
  tetgen_flags_alts = {'-Y',''};
  TT = [];
  for t = 1:numel(tetgen_flags_alts)
    tetgen_flags = tetgen_flags_alts{t};
    try
      [TV,TT,TF] = cdt(SV,SF,'TetgenFlags',tetgen_flags,'Quiet',true);
      break;
    catch
      fprintf('cdt failed with "%s"...',tetgen_flags);
    end
  end
  if isempty(TT)
    warning('Tetgen failed to maintain faces.');
    [TV,TT,TF] = cdt(SV,[],'Quiet',true);
  end
  w = winding_number(V,F,barycenter(TV,TT));
  vol = volume(TV,TT);
  avg_w = w'*(vol/sum(vol));
  if avg_w < 0
    fprintf('inside out?...');
    w = -w;
    avg_w = -avg_w;
  end
  CT = TT(w>avg_w,:);

  CF = boundary_faces(CT);           
  if surface_only
    [CV,I] = remove_unreferenced(TV,CF);
    CF = I(CF);
  end

end
