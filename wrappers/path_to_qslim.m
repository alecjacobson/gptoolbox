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
    warning([ ...
      'Dear Ladislav, is there a standard place to put executables on a pc?' ...
      'Could you put qslim there and change this accordingly?' ...
      'Thanks, Alec']);
    s = 'c:/prg/lib/qslim/Release/qslim.exe';
  elseif isunix || ismac
    % I guess this means linux
    [status, s] = system('which qslim');
    s = strtrim(s);
    if isempty(s)
      guesses = { ...
        '/usr/local/bin/qslim', ...
        '/opt/local/bin/qslim'};
      s = ...
        guesses{find(cellfun(@(guess) exist(guess,'file'),guesses),1,'first')};
      assert(~isempty(s),'Could not find qslim');
    end
  end
end
