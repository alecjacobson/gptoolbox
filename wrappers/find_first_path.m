function [s] = find_first_path(guesses,quiet)
  si = find(cellfun(@(guess) exist(guess,'file'),guesses),1,'first');
  if isempty(si)
    if nargin>1 && quiet
      s = [];
      return;
    else
      error('Could not find path');
    end
  end
  s = guesses{si};
end
