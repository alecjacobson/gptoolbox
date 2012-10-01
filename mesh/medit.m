function [s,r] = medit(V,T,F,wait)
  % MEDIT wrapper to launch medit to visualize a given tetmesh
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   T  #T by 4|5 list of tetrahedron findices (additional column is color
  %     index)
  %   F  #F by 3|4 list of face indices (additional column is color index)
  %   Optional:
  %     wait  boolean whether to wait for completion
  % Outputs:
  %   s,r  result of system call
  %

  if ~exist('wait','var')
    wait = true;
  end

  % Change these paths accordingly
  MEDIT_PATH = '/opt/local/bin/medit';
  TEMP_MESH_FILE  = '/var/tmp/temp.mesh';
  TEMP_MEDIT_FILE = '/var/tmp/temp.medit';

  % write temporary mesh
  writeMESH(TEMP_MESH_FILE,V,T,F);

  % write default medit options
  f = fopen(TEMP_MEDIT_FILE,'w');
  fprintf(f, [...
    'BackgroundColor 1 1 1\n' ...
    'LineColor 0 0 0\n'       ...
    'WindowSize 1024 800\n']);
  if size(F,2) > 3 || size(T,2) > 4
    fprintf(f,['NbMaterials\n']);
    %% collect indices
    %I = [F(:,3:end);T(:,4:end)];
    %count = max(I) - min(I);
    %fprintf(f,'%d\n',count);
    fprintf(f,'RenderMode colorshadinglines\n');
  else
    fprintf(f,'RenderMode shading + lines\n');
  end
  fclose(f);

  command = [MEDIT_PATH ' ' TEMP_MESH_FILE ' ' TEMP_MEDIT_FILE];
  if ~wait
    command = [command ' &'];
  end
  [s,r] = system(command);
end
