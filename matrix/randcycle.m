function [C,u,IA,IC] = randcycle(A)
  % RANDCYCLE Randomly reindexs a list a indices
  %
  % Inputs:
  %   A  #A list of indices with unique values u
  % Outputs:
  %   C  #C list of indices with same unique values u
  %   u,IA,IC  output of unique(A);
  %

  [u,IA,IC] = unique(A);
  % scramble u
  u = u(randperm(end));
  C = u(IC);
end
