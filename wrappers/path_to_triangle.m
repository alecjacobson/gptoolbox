function s = path_to_triangle()
  % PATH_TO_TRIANGLE Returns absolute, system-dependent path to triangle
  % executable
  %
  % Outputs:
  %   s  path to triangle as string
  %  
  % See also: triangle

  if ispc
    warning([ ...
      'Dear Ladislav, is there a standard place to put executables on a pc?' ...
      'Could you put triangle there and change this accordingly?' ...
      'Thanks, Alec']);
    s = 'c:/prg/lib/triangle/Release/triangle.exe';
  elseif isunix || ismac
    % I guess this means linux
    [status, s] = system('which triangle');
    s = strtrim(s);
    if isempty(s)
      guesses = { ...
        '/usr/local/bin/triangle', ...
        '/opt/local/bin/triangle'};
    found = find(cellfun(@(guess) exist(guess,'file'),guesses),1,'first');
    if found
      s = ...
        guesses{find(cellfun(@(guess) exist(guess,'file'),guesses),1,'first')};
    end
      assert(~isempty(s),'Could not find triangle');
    end
  end
end

