function [U,G] = remesh_botsch(V,F,varargin)
  % REMESH_BOTSCH Remesh a triangle mesh to have a desired edge length, by
  % calling gpytoolbox.remesh_botsch through matlab's python interface.
  %
  % Uses the algorithm of "A Remeshing Approach to Multiresolution Modeling"
  % [Botsch & Kobbelt 2004]: alternating iterations of subdivision, collapse,
  % edge flips and tangential smoothing.
  %
  % [U,G] = remesh_botsch(V,F)
  % [U,G] = remesh_botsch(V,F,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle mesh indices into rows of V
  %     Optional:
  %       'Iters'  followed by number of remeshing iterations {10}
  %       'EdgeLength'  followed by desired edge length {mean input edge
  %         length}
  %       'Project'  followed by whether to reproject onto the input surface
  %         (otherwise the mesh smooths over iterations) {true}
  %       'Feature'  followed by #feat list of indices into rows of V of
  %         vertices which should not move (they appear first in U, in order).
  %         Boundary vertices are always treated as features by gpytoolbox. {[]}
  % Outputs:
  %   U  #U by 3 list of output mesh vertex positions
  %   G  #G by 3 list of output triangle indices into rows of U
  %
  % Example:
  %   [V,F] = subdivided_sphere(2);
  %   [U,G] = remesh_botsch(V,F,'EdgeLength',0.1,'Iters',20);
  %   tsurf(G,U);
  %
  % See also: pip_install, pysparse, remesh_planar_patches
  %

  iters = 10;
  h = [];
  project = true;
  feature = [];
  params_to_variables = containers.Map( ...
    {'Iters','EdgeLength','Project','Feature'}, ...
    {'iters','h','project','feature'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  assert(size(V,2)==3,'Only 3d vertex positions supported');
  assert(size(F,2)==3,'Only triangle meshes supported');

  % Matlab passes a bare matrix to python as a memoryview and a bare vector as
  % an array.array; py.numpy.array turns either into a correctly shaped
  % ndarray. atleast_2d/atleast_1d undo the squeezing of single-row inputs. The
  % classes matter: the c++ binding demands float64 positions and int32
  % 0-indexed indices.
  pyV = py.numpy.atleast_2d(py.numpy.array(double(V)));
  pyF = py.numpy.atleast_2d(py.numpy.array(int32(F)-1));
  pyfeature = py.numpy.atleast_1d(py.numpy.array(int32(feature(:)')-1));
  if isempty(h)
    % None tells gpytoolbox to use the mean input edge length
    pyh = py.None;
  else
    pyh = double(h);
  end

  try
    out = py.gpytoolbox.remesh_botsch(pyV,pyF,int32(iters),pyh,logical(project),pyfeature);
  catch ME
    assert_gpytoolbox();
    rethrow(ME);
  end

  U = double(out{1});
  G = double(out{2})+1;
end

function assert_gpytoolbox()
  % Issue an actionable error if gpytoolbox is missing from the python
  % environment matlab is currently using (otherwise return silently).
  try
    py.importlib.import_module('gpytoolbox');
    return;
  catch
  end
  pe = pyenv;
  if strlength(pe.Version) == 0
    error('remesh_botsch:no_python', ...
      ['Matlab has no python environment configured.\n' ...
       'Point it at an interpreter with gpytoolbox installed, e.g.:\n\n' ...
       '    pyenv(''Version'',''/usr/local/bin/python3'')\n']);
  end
  error('remesh_botsch:no_gpytoolbox', ...
    ['gpytoolbox is not installed in the python environment matlab is using\n' ...
     '(%s).\nInstall it there with:\n\n' ...
     '    pip_install(''gpytoolbox'')\n\n' ...
     'or switch interpreters with `pyenv(''Version'',''/path/to/python3'')`.\n'], ...
    pe.Executable);
end
