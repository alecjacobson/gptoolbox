function [Y,I] = maxnz(X)
  % MAXNZ  find maximum nonzero entry in columns of X
  %
  % [Y,I] = maxnz(X)
  %
  % Inputs:
  %   X  m by n sparse matrix
  % Outputs
  %   Y  n list of maximum non zero entries in each column of X
  %   I  n list of row indices to maximum non-zero entries in each column of X
  %
  % See also: min, max, minnz
  %

  [XI,XJ,XV] = find(X);
  X_i = sparse(XI,XJ,XV.^-1,size(X,1),size(X,2));
  [minX_i,minI_i] = min(X_i);
  [maxX,maxI] = max(X);
  minX_i(minX_i==0) = -inf;
  maxX(maxX==0) = -inf;
  Y = [minX_i.^-1;maxX];
  [Y,J] = max(Y);
  I = [minI_i;maxI];
  I = I(sub2ind(size(I),J,1:size(I,2)));
end

