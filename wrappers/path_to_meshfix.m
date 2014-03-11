function s = path_to_meshfix()
  % PATH_TO_MESHFIX Returns absolute, system-dependent path to meshfix executable
  %
  % Outputs:
  %   s  path to meshfix as string
  %  
  % See also: meshfix

  if ispc
    warning([ ...
      'Dear Ladislav, is there a standard place to put executables on a pc?' ...
      'Could you put meshfix there and change this accordingly?' ...
      'Thanks, Alec']);
    s = 'c:/prg/lib/meshfix/Release/meshfix.exe';
  elseif ismac
    s = '/usr/local/igl/libigl/external/MeshFix/meshfix';
  elseif isunix 
    % I guess this means linux
    s = '/usr/local/bin/meshfix';
  end
end
