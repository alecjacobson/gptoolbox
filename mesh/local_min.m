function [min_v,min_i,min_diff] = local_min(F,S,ep)
  % LOCAL_MIN  find values and indices of local minima of a scalar field
  % defined on a mesh
  % 
  % [min_v,min_i,min_diff] = local_min(F,S,ep)
  %
  % Inputs:
  %   F #faces by 3 list of triangle indices
  %   S #vertices by 1 list of scalar values
  % Outputs:
  %   min_v  #minima list of values at local minima
  %   min_i  #minima list of indices into S of local minima
  %   min_diff #minima list of min difference between each minima and its
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
  minA = minnz(A)';
  [min_i,~] = find(minA>-ep);
  min_diff = minA(min_i);
  min_v = S(min_i);
    
end

