function s = path_to_meshfix()
  % PATH_TO_MESHFIX Returns absolute, system-dependent path to meshfix executable
  %
  % Outputs:
  %   s  path to meshfix as string
  %  
  % See also: meshfix

  if ispc
    s = 'c:/prg/lib/meshfix/Release/meshfix.exe';
  elseif isunix || ismac
    % I guess this means linux
    [status, s] = system('which meshfix');
    s = strtrim(s);
    if status ~= 0
      guesses = { ...
        '/usr/local/bin/meshfix', ...
        '/opt/local/bin/meshfix', ...
        '/usr/local/igl/libigl/external/MeshFix/meshfix', ...
        '/usr/local/libigl/external/MeshFix/meshfix'};
      s = find_first_path(guesses);
    end
  end
end
