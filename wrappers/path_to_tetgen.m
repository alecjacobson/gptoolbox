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
    [~,s] = system('which tetgen');
    s = strtrim(s);
    if isempty(s)
      guesses = { ...
        '/usr/local/bin/tetgen', ...
        '/opt/local/bin/tetgen', ...
        '/usr/local/igl/libigl/external/tetgen/tetgen', ...
        [path_to_libigl '/external/tetgen/tetgen'], ...
        '/usr/local/libigl/external/tetgen/tetgen'};
      s = ...
        guesses{find(cellfun(@(guess) exist(guess,'file'),guesses),1,'first')};
      assert(~isempty(s),'Could not find tetgen');
    end
  end
end
