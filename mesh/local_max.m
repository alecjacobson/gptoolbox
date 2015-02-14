function [max_v,max_i,max_diff] = local_max(F,S,ep)
  % LOCAL_MAX  find values and indices of local maxima of a scalar field
  % defined on a mesh
  %
  % [max_v,max_i,max_diff] = local_max(F,S,ep)
  %
  % Inputs:
  %   F #faces by 3 list of triangle indices
  %   S #vertices by 1 list of scalar values
  % Outputs:
  %   max_v  #maxima list of values at local maxima
  %   max_i  #maxima list of indices into S of local maxima
  %   max_diff #maxima list of max difference between each maxima and its
  %     neighbors
  %

  if ~exist('ep','var')
    ep = 0;
  end

  % number of vertices
  n = numel(S);
  if size(F,2) == 2
    E = F;
  else
    E = edges(F);
  end
  E = [E;fliplr(E)];
  A = sparse( ...
    E(:,1), ...
    E(:,2), ...
    S(E(:,1))-S(E(:,2)), ...
    n,n);
  maxA = maxnz(A)';
  [max_i,~] = find(maxA<ep);
  max_diff = maxA(max_i);
  max_v = S(max_i);
    
end
