function [W,G] = meshfix(V,F) 
  % Call Marco Atene's MESHFIX on a mesh
  %
  % [W,G] = meshfix(V,F) 
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices
  % Outputs:
  %   W  #W by 3 list of output mesh vertex positions
  %   G  #G by 3 list of output triangle indices
  %
  prefix = tempname;
  off_filename = [prefix '.off'];
  off_filename_fixed = [prefix '_fixed.off'];
  writeOFF(off_filename,V,F);
  flags = '';
  command = [path_to_meshfix ' ' flags ' ' off_filename];
  [status, result] = system(command);
  if status ~= 0
    fprintf(command);
    error(result);
  end
  [W,G] = readOFF(off_filename_fixed);

  delete(off_filename);
  delete(off_filename_fixed);
end
