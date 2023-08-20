function s = path_to_triangle()
  % PATH_TO_TRIANGLE Returns absolute, system-dependent path to triangle
  % executable
  %
  % Outputs:
  %   s  path to triangle as string
  %  
  % See also: triangle

  if ispc
    s = 'c:/prg/lib/triangle/Release/triangle.exe';
  elseif isunix || ismac
    [status, s] = system('which triangle');
    s = strtrim(s);
    if status ~= 0
      guesses = { ...
        '/usr/local/bin/triangle', ...
        '/opt/local/bin/triangle'};
      s = find_first_path(guesses);
    end
  end
end

