function y = cumoccurrence(x)
  % y = cumoccurrence(x)
  %
  % Inputs:
  %   x  list of elements
  % Outputs
  %   y  list of cummulative occurances, so that y(i) = j indicates that x(i)
  %     is the jth occurance of the value x(i) so far.
  %

  % Slightly faster than mine. Really the index-tracking sort is dominating.
  % https://stackoverflow.com/a/64034583/148668
  [s, is] = sort(x);
  d = [1 ;diff(s)];
  f = find(d);
  d(f) = f;
  ic = cummax(d);
  y = zeros(size(x));
  y(is) = (2 : numel(s) + 1).' - ic;
end
