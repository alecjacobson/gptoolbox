function s = path_to_convert()
  % PATH_TO_CONVERT
  %
  % s = path_to_convert()
  %
  % Outputs:
  %   s path to convert executable
  %
  % See also: path_to_qslim
  %

  if ispc
    % replace this with path
    s = 'c:/prg/lib/convert/Release/convert.exe';
  elseif isunix || ismac
    [status, s] = system('which convert');
    s = strtrim(s);
    if status ~= 0
      guesses = { ...
        '/usr/local/bin/convert', ...
        '/opt/local/bin/convert'};
      s = find_first_path(guesses);
    end
  end
end

