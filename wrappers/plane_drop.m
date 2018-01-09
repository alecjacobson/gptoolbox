function [x,R,v,omega,M,I,dt,t] = ...
  plane_drop(VV,FF,geo_idx,x0,R0,v0,omega0,rho,fixed,varargin)
  % PLANE_DROP  Drop some shit on the ground using SCISIM.
  %
  % This assumes you've installed scisim and implicittoolkit. Congratulations if
  % you've installed scisim and implicittoolkit!
  %
  % [x,R,v,omega,M,I,dt,t] = plane_drop(VV,FF,geo_idx,x,R,v,omega,rho,fixed)
  %
  % Inputs:
  %   VV  ngeo list of #VV{i} by 3 lists of mesh vertex positions
  %   FF  ngeo list of #FF{i} by 3 list of triangle mesh indices into VV{i}
  %   geo_idx  k list of indices into VV,FF indicate which mesh to use for each
  %     rigid body (first index is 1 like matlab)
  %   x0  k by 3 list of initial center of mass locations
  %   R0  k by 3 list of initial rotations (Euler angles), [] --> zero
  %   v0  k by 3 list of initial velocities, [] --> zero
  %   omega0  k by 3 list of initial angular velocities (Euler angles?), 
  %     [] --> zero
  %   rho  k list of densities(?) [] --> one
  %   fixed  k list of flags whether object is fixed [] --> false
  %   Optional:
  %     'CellWidth'  followed by scalar or ngeo list of cell width value(s) to
  %       use for signed distance field {0.1*bbd{i}}
  %     'Duration'  followed by simulation duration in seconds {1}
  %     'FPS'  followed by frames-per second rate of output animation {30}
  %     'GridPadding'  followed by scalar or ngeo list of grid padding value(s) to
  %       use for signed distance field 
  % Outputs
  %   x  k by 3 by #frames list of center of mass locations
  %   R  k by 9 by #frames list of rotations (matrices...)
  %   v  k by 3 by #frames list of velocities
  %   omega  k by 3 by #frames list of angular velocities (Euler angles?)
  %   M  k by 9 by #frames list of moments?
  %   I  k by 9 by #frames list of moments?
  %   dt  #frames list of timesteps (should be constant)
  %   t  #frames list of times (should be equal to cumsum(dt))
  %   
  %   

  function thisElement = append_with_attributes(docNode,name,attr)
    thisElement = docNode.createElement(name);
    for i = 1:2:numel(attr)
      if ischar(attr{i+1})
        str = attr{i+1};
      else
        str = strtrim(sprintf('%g ',attr{i+1}'));
      end
      thisElement.setAttribute(attr{i},str);
    end
    docRootNode.appendChild(thisElement);
  end

  % default values
  cell_width = [];
  end_time = 1;
  fps = 30;
  grid_padding = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'CellWidth','Duration','FPS','GridPadding'}, ...
    {'cell_width','end_time','fps','grid_padding'});
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

  if ~iscell(VV)
    VV = {VV};
  end
  if ~iscell(FF)
    FF = {FF};
  end
  ngeo = numel(VV);
  if isempty(cell_width)
    cell_width = cellfun(@(V) 0.1*normrow(max(V)-min(V)),VV);
  end
  if isempty(grid_padding)
    grid_padding = cellfun(@(V) 0.1*normrow(max(V)-min(V)),VV);
  end

  % generate h5 file for each input mesh
  h5 = cell(ngeo,1);
  for i = 1:ngeo
    path_to_implicittoolkit_cli = '/Users/ajx/Dropbox/implicittoolkit/build/implicittoolkitcli/implicittoolkit_cli';
    obj = sprintf('plane_drop-%d.obj',i-1);
    h5{i} = sprintf('%s/plane_drop-%d.h5',pwd,i-1);
    writeOBJ(obj,VV{i},FF{i});
    cmd = sprintf( ...
      '%s -w %0.17f -p %0.17f %s %s', ...
      path_to_implicittoolkit_cli, ...
      cell_width(min(i,end)),grid_padding(min(i,end)),obj,h5{i});
    [status,r] = system(cmd);
    if status ~= 0
      error(sprintf('%s\n\n%s',cmd,r));
    end
    delete(obj);
  end

  k = numel(geo_idx);
  assert(min(geo_idx)>=1 && max(geo_idx)<=ngeo);
  if isempty(R0)
    R0 = zeros(k,3);
  end
  if isempty(v0)
    v0 = zeros(k,3);
  end
  if isempty(omega0)
    omega0 = zeros(k,3);
  end
  if isempty(rho)
    rho = ones(k,1);
  end
  if isempty(fixed)
    fixed = false(k,1);
  end

  docNode = com.mathworks.xml.XMLUtils.createDocument('rigidbody3d_scene');
  docRootNode = docNode.getDocumentElement;
  append_with_attributes(docNode,'camera_perspective',{ ...
    'theta','1.3172', ...
    'phi','0.265398', ...
    'rho','32.6039', ...
    'lookat','5.83134 5.26439 1.44123', ...
    'up','0 1 0', ...
    'fps','60', ...
    'render_at_fps','1', ...
    'locked','0'});

  append_with_attributes(docNode,'integrator',{'type','dmv','dt','1/10800'});

  append_with_attributes(docNode,'sobogus_friction_solver',{ ...
    'mu','0.4', ...
    'CoR','0.3', ...
    'max_iters','5000', ...
    'tol','5.0e-8', ...
    'eval_every','50', ...
    'staggering','geometric'});

  add_gravity = true;
  if add_gravity
    append_with_attributes(docNode,'near_earth_gravity',{'f','0 -981 0'});
  end

  add_plane = true;
  plane_center = [0.0 0.0 0.0];
  plane_normal = [0.0 1.0 0.0];
  if add_plane
    append_with_attributes(docNode,'static_plane', ...
      {'x',plane_center,'n',plane_normal});
    append_with_attributes(docNode,'static_plane_renderer', ...
      {'plane','0','r',[20 20]});
  end

  for i = 1:ngeo
    append_with_attributes(docNode,'geometry',{'type','mesh','filename',h5{i}});
  end


  for i=1:k
    append_with_attributes(docNode,'rigid_body_with_density',{ ... 
      'x',x0(i,:), ...
      'R',R0(i,:), ...
      'v',v0(i,:), ...
      'omega',omega0(i,:), ...
      'rho',rho(i), ...
      'fixed',fixed(i), ...
      'geo_idx',geo_idx(i)-1});
  end
  xmlFileName = 'plane_drop.xml';
  xmlwrite(xmlFileName,docNode);
  %type(xmlFileName);
  path_to_rigidbody3d_cli = ...
    '/Users/ajx/Dropbox/scisim/build/rigidbody3dcli/rigidbody3d_cli';
  path_to_rigidbody3d_qt4 = ...
    '/Users/ajx/Dropbox/scisim/build/rigidbody3dqt4/rigidbody3d_qt4.app/Contents/MacOS/rigidbody3d_qt4';

  if ~exist('plane_drop-data','dir')
    mkdir('plane_drop-data/');
  else
    % what a dangerous command to have lying in here
    delete('plane_drop-data/*.h5');
  end

  % while you wait, why don't you preview it in qt?
  qt = sprintf('%s %s &',path_to_rigidbody3d_qt4,xmlFileName);
  %fprintf(qt);
  [status,r] = system(qt);
  cmd = sprintf( ...
    '%s -o %s -e %0.17f -f %0.17f %s', ...
    path_to_rigidbody3d_cli,'plane_drop-data/',end_time,fps,xmlFileName);
  [status,r] = system(cmd);
  if status ~= 0
    error(sprintf('\n%s\n\n%s',cmd,r));
  else
    % Erase qt line
  %  fprintf(repmat('\b',[1 length(qt)]));
  end
  


  h5s = dir('plane_drop-data/*.h5');
  m = numel(h5s);
  x     = zeros(k,3,m);
  R     = zeros(k,9,m);
  v     = zeros(k,3,m);
  omega = zeros(k,3,m);
  M     = zeros(k,1,m);
  I     = zeros(k,9,m);
  fixed = zeros(k,1,m);
  dt = zeros(m,1);
  t = zeros(m,1);
  for hi = 1:numel(h5s)
    [iter,dt(hi),t(hi), ...
      x(:,:,hi), ...
      R(:,:,hi), ...
      v(:,:,hi), ...
      omega(:,:,hi), ...
      M(:,:,hi), ...
      I(:,:,hi)] = ...
      readSCISIM([h5s(hi).folder '/' h5s(hi).name]);
  end


end
