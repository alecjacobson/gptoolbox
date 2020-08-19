function enable_toolbox(name)
  local_toolbox_dir = ['./' name];
  if exist(local_toolbox_dir,'dir');
     rmdir(local_toolbox_dir,'s');
  end
  mkdir(local_toolbox_dir);
  mkdir([local_toolbox_dir '/private'])
  addpath(local_toolbox_dir) % important to do this before copying
  matlabroot = '/Applications/MATLAB_R2018a.app/';
  copyfile( ...
    [matlabroot '/toolbox/' name '/' name '/*.m'], ...
    [local_toolbox_dir '/'])
  copyfile( ...
    [matlabroot '/toolbox/' name '/' name '/private/*.m'], ...
    [local_toolbox_dir '/private'])
end
