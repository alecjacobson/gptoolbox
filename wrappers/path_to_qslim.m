function s = path_to_qslim()
  % PATH_TO_QSLIM
  %
  % s = path_to_qslim()
  %
  % Outputs:
  %   s path to qslim executable
  %
  % See also: qslim
  %

  if ispc
    s = 'c:/prg/lib/qslim/Release/qslim.exe';
  elseif isunix || ismac
    % I guess this means linux
    [status, s] = system('which qslim');
    s = strtrim(s);
    if status ~= 0
      guesses = { ...
        '/usr/local/bin/qslim', ...
        '/opt/local/bin/qslim'};
      s = find_first_path(guesses);
    end
  end
end
