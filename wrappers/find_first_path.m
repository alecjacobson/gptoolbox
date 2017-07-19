function [s] = find_first_path(guesses)
  si = find(cellfun(@(guess) exist(guess,'file'),guesses),1,'first');
  assert(~isempty(si),'Could not find path');
  s = guesses{si};
end
