function U = composite_harmonic_mapping(V,F,b,bc,varargin)
  % COMPOSITE_HARMONIC_MAPPING Compute a "composite harmonic mapping in the
  % plane. This is a similar in spirit to "(bijective) composite mean value
  % mappings" in [Schneider et al.] and also (iff outer_iters=1 and bc=V(b,:)
  % then this is equivalent to "Embedding a triangular graph within a given
  % boundary" [Xu et al. 2011].
  %
  % U = composite_harmonic_mapping(V,F,b,bc,varargin)
  %
  % Inputs:
  %   V  #V by 2 list of vertices
  %   F  #F by 2 list of triangle indices
  %   b  #b list of indices of V of constrained vertices
  %   bc  #bc by 2 list of boundary conditions
  %   Optional:
  %     'OuterIters'  followed by number of iterations for the morphing {20}
  %     'InnerIters'  followed by number of iterations per morph-iteration {2}
  % Outputs:
  %   U  #V by 2 list of output vertex positions
  %
  % Example:
  %   [V,F,VT,FT] = readOBJ('~/Dropbox/models/textured-models-from-daniele/animal/animal.obj');
  %   [UE,UT,OT] = straighten_seams(V,F,VT,FT,'Tol',inf);
  %   b = unique(OT);
  %   bc = UT(b,:);
  %   WCH = composite_harmonic_mapping(VT,FT,b,bc);

  outer_iters = 20;
  inner_iters = 2;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'OuterIters','InnerIters'}, ...
    {'outer_iters','inner_iters'});
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

  % TODO: This should be replaced with something better like the 2015 Lipman
  % curve morphing thing.
  lerp_bc = @(t) V(b,:)+t.*(bc-V(b,:));
  func_bc = lerp_bc;
  U = V;
  for outer_iter = 1:outer_iters
    t = outer_iter/outer_iters;
    for inner_iter = 1:inner_iters
      % TODO: if the delta between each iteration is small then this might
      % benefit from switching to an iterative solver using a warm start: pcg
      U = kharmonic(U,F,b,func_bc(t),1);
    end
  end

end
