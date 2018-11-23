function combs = nmultichoosek(values, k)
  % NMULTICHOOSEK Generate all k-combinations with repetitions of a given
  % set.
  %
  % Inputs:
  %   values  #values set to consider (should already contain unique elements)
  %   k  size of combination to draw
  % Outputs:
  %   combs  nchoosek(#values,k) by k list of combinations
  %
  % https://stackoverflow.com/a/28284672/148668
  assert(numel(values)>=k);
  n = numel(values);
  combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
  combs = reshape(values(combs),[],k);
end
