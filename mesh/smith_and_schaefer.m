function [U,data] = smith_and_schaefer(V,F,U0,varargin)
  % SMITH_AND_SCHAEFER
  %
  % U = smith_and_schaefer(V,F,U0,varargin)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indices into rows of V
  %   U0  #V by 2 list of initial feasible state positions
  %   Optional:
  %     'Data' see output
  %     'MaxIters' followed by maximimum number of iterations {100}
  %     'Quiet' followed by whether to not print convergence info
  %     'Tol'  followed by boundary barrier tolerance  see
  %       self_collision_barrier.m
  % Outputs:
  %   U  #V by 2 list of output feasible state positions
  %   data  struct containing information so that restarting iterations is
  %     continuous.
  %

  max_iters = 100;
  tol = [];
  data = struct();
  quiet = false;
  prevent_boundary_overlaps = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'PreventBoundaryOverlaps','MaxIters','Quiet','Tol','Data'}, {'prevent_boundary_overlaps','max_iters','quiet','tol','data'});
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

  % use tutte embedding if no initial feasible state
  if isempty(U0)
    U0 = tutte(V,F);
  end

  % default boundary barrier tolerance
  O = boundary_faces(F);
  if isempty(tol)
    tol = avgedge(V,O)/4;
  end

  if ~isfield(data,'u') || isempty(data.u)
    data.u = 1e-1;
  end
  if ~isfield(data,'iter') || isempty(data.iter)
    data.iter = 0;
  end
  if ~isfield(data,'iters_since_stall') || isempty(data.iters_since_stall)
    data.iters_since_stall = 0;
  end
  if ~isfield(data,'K') || isempty(data.K)
    data.K = 10;
  end

  U = U0;
  % ∞ if intersecting, 0 otherwise
  self_intersecting = ...
   @(U) 1/(~is_self_intersecting(reshape(U,[],2),O))-1;
  % ∞ if flipped, 0 otherwise
  flipped = @(X)  max((1./(doublearea(reshape(X,[],2),F)>0) - 1));
  assert(~prevent_boundary_overlaps || self_intersecting(U) == 0);
  assert(flipped(U) == 0);
  % u can be updated multiplicatively /=2  *=2
  u2w = @(u) u/(u+1);

  data.converged = false;
  for iter = data.iter+(1:max_iters)
    data.iter = iter;
    U0 = reshape(U,[],2);
    [fsd,Gsd,Hsd] = symmetric_dirichlet(U0,F,V);
    if prevent_boundary_overlaps
      [fb,Gb,Hb] = self_collision_barrier(U0,O,tol);
      f0 = fsd + fb;
      G0 = Gsd + Gb;
      H0 = Hsd + Hb;
    else
      f0 = fsd;
      G0 = Gsd;
      H0 = Hsd;
    end
    while true
      stalled = false;
      w = u2w(data.u);
      Htilde = ((1-w)*H0+w*speye(size(H0)));
      %Dec = decomposition(Htilde,'chol','upper');
      %dU = Dec\(-G0);
      dU = Htilde\(-G0);
  
      % Pure gradient descent
      %[fsd,Gsd] = symmetric_dirichlet(U0,F,V);
      %[fb,Gb] = self_collision_barrier(U0,O,tol);
      %f0 = fsd + fb;
      %G0 = Gsd + Gb;
      %dU = -G0;

      if max(abs(dU))>1e-10
        if prevent_boundary_overlaps
          f_funs = @(X) [
            symmetric_dirichlet(   reshape(X,[],2),F,V) ...
            self_collision_barrier(reshape(X,[],2),O,tol) ...
            flipped(X) ...
            self_intersecting(X)];
        else
          f_funs = @(X) [
            symmetric_dirichlet(   reshape(X,[],2),F,V) ...
            flipped(X) ];
        end
        f_fun = @(X) sum(f_funs(X));
        [t,U1,data.f] = backtracking_line_search(f_fun,U0(:),G0,dU,0.01,0.5);
        if f0-data.f < 1e-15
          if ~quiet
            fprintf('    small objective change\n');
          end
          stalled = true;
        end
        if ~stalled && t < 1e-10
          if ~quiet
            fprintf('    small step size\n');
          end
          stalled = true;
        end
      else
        if ~quiet
          fprintf('    small step direction\n');
        end
        stalled = true;
      end
      if stalled 
        if (1-u2w(data.u))<1e-10
          if ~quiet
            fprintf('Stalled with gradient descent. Game over.\n');
          end
          data.converged = true;
          break;
        end
        data.u = data.u*2;
        data.iters_since_stall = 0;
        if ~quiet
          fprintf('    u: %g\n',data.u);
        end
        continue;
      end
      data.iters_since_stall = data.iters_since_stall+1;
      U = reshape(U1,[],2);
      if ~quiet
        fprintf('% 4d.%02d: %g\n',data.iter,data.iters_since_stall,data.f);
      end
      if data.iters_since_stall>= data.K
        data.u = max(data.u/2,1e-8);
        data.iters_since_stall = 0;
        if ~quiet
          fprintf('    u: %g\n',data.u);
        end
      end
      break;
    end
    if data.converged
      break;
    end
  
  end


end
