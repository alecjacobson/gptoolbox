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
    warning([ ...
      'Dear Ladislav, is there a standard place to put executables on a pc?' ...
      'Could you put convert there and change this accordingly?' ...
      'Thanks, Alec']);
    s = 'c:/prg/lib/convert/Release/convert.exe';
  elseif isunix || ismac
    % I guess this means linux
    [status, s] = system('which convert');
    s = strtrim(s);
    if isempty(s)
      guesses = { ...
        '/usr/local/bin/convert', ...
        '/opt/local/bin/convert'};
      E = cellfun(@(guess) exist(guess,'file'),guesses);
      assert(any(E),'Could not find convert');
      s = guesses{find(E,1,'first')};
    end
  end
end

