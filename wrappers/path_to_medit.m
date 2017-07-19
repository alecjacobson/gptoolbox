function [s] = path_to_medit()
  % PATH_TO_MEDIT Return path to medit executable
  %
  % s = path_to_medit()
  %
  % Outputs:
  %   s path to medit executable
  %
  % See also: medit
  %


  if ispc
    warning([ ...
      'Dear Ladislav, is there a standard place to put executables on a pc?' ...
      'Could you put medit there and change this accordingly?' ...
      'Thanks, Alec']);
    s = 'c:/prg/lib/medit/Release/medit.exe';
  elseif isunix || ismac
    % I guess this means linux
    [status, s] = system('which medit');
    s = strtrim(s);
    if isempty(s)
      guesses = { ...
        '/usr/local/bin/medit', ...
        '/opt/local/bin/medit'};
      s = find_first_path(guesses);
    end
  end
end
