function s = path_to_tetgen()
  % PATH_TO_TETGEN Returns absolute, system-dependent path to tetgen executable
  %
  % Outputs:
  %   s  path to tetgen as string
  %  
  % See also: tetgen

  if ispc
    % replace this with path
    s = 'c:/prg/lib/tetgen/Release/tetgen.exe';
  elseif ismac || isunix
    [status,s] = system('which tetgen');
    s = strtrim(s);
    if status ~= 0
      guesses = { ...
        '/usr/local/bin/tetgen', ...
        '/opt/local/bin/tetgen', ...
        '/Users/ajx/Repos/tetgen/build/tetgen', ...
        '/usr/local/igl/libigl/external/tetgen/tetgen', ...
        '/usr/local/libigl/external/tetgen/tetgen'};
      s = find_first_path(guesses);
    end
  end
end
