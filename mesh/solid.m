function [CV,CT,TV,TT,w,SV,SF] = solid(V,F,varargin)
  % SOLID Clean a mesh using the generalized winding number pipeline with some
  % default settings [Jacobson et al. 2013]: "Shoot for the moon!!!!"
  %
  % [CV,CT] = solid(V,F)
  % [CV,CT] = solid(V,F,'ParameterName',ParameterValue,..)
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
  [SV,SF] = clean(V,F,'Quiet',true,'Single',true);
  tetgen_flags_alts = {'-Y','-Y -T1e-16','','-T1e-16'};
  TT = [];
  for t = 1:numel(tetgen_flags_alts)
    tetgen_flags = tetgen_flags_alts{t};
    try
      [TV,TT] = cdt(SV,SF,'TetgenFlags',tetgen_flags,'Quiet',true);
      fprintf('cdt succeeded with "%s"...\n',tetgen_flags);
      break;
    catch
      fprintf('cdt failed with "%s"...\n',tetgen_flags);
    end
  end
  if isempty(TT)
    error('All tetgen flags alternatives failed');
  end
  w = winding_number(V,F,barycenter(TV,TT));
  vol = volume(TV,TT);
  avg_w = w'*(vol/sum(vol));
  if avg_w < 0
    fprintf('inside out?...\n');
    w = -w;
    avg_w = -avg_w;
  end
  CT = TT(w>avg_w,:);

  if surface_only
    CT = boundary_faces(CT);           
  end
  [CV,I] = remove_unreferenced(TV,CT);
  CT = I(CT);
  
end
