function [Y,I] = minnz(X)
  % MINNZ  find minimum nonzero entry in columns of X
  %
  % [Y,I] = minnz(X)
  %
  % Inputs:
  %   X  m by n sparse matrix
  % Outputs
  %   Y  n list of minimum non zero entries in each column of X
  %   I  n list of row indices to minimum non-zero entries in each column of X
  %
  % See also: min, max, maxnz
  %
  %
  
  [XI,XJ,XV] = find(X);
  X_i = sparse(XI,XJ,XV.^-1,size(X,1),size(X,2));
  [maxX_i,maxI_i] = max(X_i);
  [minX,minI] = min(X);
  maxX_i(maxX_i==0) = inf;
  minX(minX==0) = inf;
  Y = [maxX_i.^-1;minX];
  [Y,J] = min(Y);
  I = [maxI_i;minI];
  I = I(sub2ind(size(I),J,1:size(I,2)));
end
