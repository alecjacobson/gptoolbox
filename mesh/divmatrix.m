function D = divmatrix(V,F)
% D = gradmatrix(V,F)
% V:  #V x 3 matrix of vertex coordinates
% F:  3 x #F  matrix of indices of triangle corners
% returns  #V x #V matrix of grad weights
  warning('Obsolete, use `div` instead');

  FT = F';
  E = ...
    unique( ...
      sort([FT(:,1) FT(:,2); FT(:,2) FT(:,3); FT(:,3) FT(:,1)]')', ...
      'rows');
  norms = sqrt(sum((V(E(:,1),:) - V(E(:,2),:)).^2,2));  
  i = [E(:,1) E(:,2) E(:,1) E(:,2)];
  j = [E(:,2) E(:,1) E(:,1) E(:,2)];
  v = [norms norms -norms -norms] ;
  % for repeated indices (i,j) sparse automatically sums up elements, as we want
  D = sparse(i,j,v,size(V,1),size(V,1));
end
