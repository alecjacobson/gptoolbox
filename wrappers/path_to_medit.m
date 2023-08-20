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
    s = 'c:/prg/lib/medit/Release/medit.exe';
  elseif isunix || ismac
    % I guess this means linux
    [status, s] = system('which medit');
    s = strtrim(s);
    if status ~= 0
      guesses = { ...
        '/usr/local/bin/medit', ...
        '/opt/local/bin/medit'};
      s = find_first_path(guesses);
    end
  end
end
