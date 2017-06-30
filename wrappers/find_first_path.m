function s = find_first_path(guesses)
      s = ...
        guesses{find(cellfun(@(guess) exist(guess,'file'),guesses),1,'first')};
      assert(~isempty(s),'Could not find path');
end
