function s = path_to_tetgen()
  % PATH_TO_TETGEN Returns absolute, system-dependent path to tetgen executable
  %
  % Outputs:
  %   s  path to tetgen as string
  %  
  % See also: tetgen

  if ispc
    warning([ ...
      'Dear Ladislav, is there a standard place to put executables on a pc?' ...
      'Could you put tetgen there and change this accordingly?' ...
      'Thanks, Alec']);
    s = 'c:/prg/lib/tetgen/Release/tetgen.exe';
  elseif ismac
    s = '/usr/local/bin/tetgen';
  elseif isunix 
    % I guess this means linux
    s = '/usr/local/bin/tetgen';
  end
end
