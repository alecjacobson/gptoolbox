function U = composite_mvc(V,P0,E,P1,varargin)
  % COMPOSITE_MVC A simplified implementation of "Bijective Composite Mean
  % Value Mappings" in [Schneider et al.]. Here, the polygon morph is a linear
  % interpolation with uniform temporal steps.
  %
  % U = composite_mvc(V,P0,E,P1,varargin)
  %
  % Inputs:
  %   V  #V by 2 list of points to be tracked
  %   P0  #P by 2 list of corner positions at t=0
  %   E  #E by 2 list of oriented edge indices into P0/P1
  %   P1  #P by 2 list of corner positions at t=1
  %   Optional:
  %     'OuterIters'  followed by number of iterations for the morphing {20}
  % Outputs:
  %   U  #V by 2 list of output vertex positions
  %
  % Example:
  %   [V,F,VT,FT] = readOBJ('~/Dropbox/models/textured-models-from-daniele/animal/animal.obj');
  %   [UE,UT,OT] = straighten_seams(V,F,VT,FT,'Tol',inf);
  %   [P0,J,I] = remove_unreferenced(VT,OT);
  %   P1 = UT(I,:);
  %   E = J(OT);
  %   WCMV = composite_mvc(VT,P0,E,P1);
  %   WCMV(I,:) = P1;
  %

  outer_iters = 20;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'OuterIters'}, ...
    {'outer_iters'});
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
  lerp_bc = @(t) P0+t.*(P1-P0);
  func_bc = lerp_bc;

  U = V;
  bc_t = func_bc(0);
  for outer_iter = 1:outer_iters
    t = outer_iter/outer_iters;
    W = mvc(U,bc_t,E,'EnforceBoundaryConditions',false);
    % Update control polygon positions
    bc_t = func_bc(t);
    % Map query points
    U = W*bc_t;
  end
end
