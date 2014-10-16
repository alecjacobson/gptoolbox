function cspy(S,epsilon)
  % cspy Visualize sparsity pattern. coloring positive and negative entries of S
  % in red and blue respectively
  %
  % cspy(S)
  % cspy(S,epsilon)
  %
  % Inputs:
  %   S  m by n (sparse) matrix
  %   epsilon  zero value {0}
  %
  % See also: spy
  %

  if nargin < 2
    epsilon = 0;
  end
  spy(S>epsilon,'r');
  hold on;
  spy(S<-epsilon,'b');
  hold off;
end
