function clean_dir(dir_arg)
% CLEAN_DIR Clean an entire directory of models
%
%
% Inputs:
%   dir_arg  argument for `dir` like '*.OBJ'
%

%addpath ../
files = dir(dir_arg);
for f = 1:numel(files);
  file = files(f).name;
  [path,name,ext] = fileparts(file);
  fprintf('%s...',file);
  [V,F] = load_mesh(file,'Quiet',true);
  mesh_fix_success = true;
  try
    [CV,CF] = meshfix(V,F);
  catch
    fprintf('meshfix crashed...');
    mesh_fix_success = false;
  end
  if size(CF,1) < 0.5*size(F,1)
    fprintf('meshfix removed more than half...');
    tsurf(CF,CV);
    axis equal;
    title(file,'FontSize',15);
    drawnow;
    mesh_fix_success = false;
  end
  if mesh_fix_success
    output_name = fullfile(path,[name '_fixed.obj']);
  else
    try
      [CV,CF] = winding_number_clean(V,F,'SurfaceOnly',true);
      output_name = fullfile(path,[name '_solid.obj']);
    catch(e)
      fprintf('failure.\n');
      continue;
    end
  end
  try
    [CCV,CCT,CCF] = tetgen(CV,CF);
    fprintf('success!\n');
  catch
    fprintf('tetgen failed on CV,CF...');
    fprintf('semi-success?\n');
  end
  tsurf(CF,CV);
  axis equal;
  title(file,'FontSize',15);
  drawnow;
  writeOBJ(output_name,CV,CF);
end

end
