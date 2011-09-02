function G = gradmatrix(V,F,N)
% G = gradmatrix(V,F)
% V:  #V x 3 matrix of vertex coordinates
% F:  3 x #F  matrix of indices of triangle corners
% N:  #N x 1 list of indices for which we want G
% returns  #E x #V matrix of grad equations
% where #E is the number of edges in region of N

  FT = limit_faces(F',N,1);
  E = ...
    unique( ...
      sort([FT(:,1) FT(:,2); FT(:,2) FT(:,3); FT(:,3) FT(:,1)]')', ...
      'rows');
  i = [1:size(E,1) 1:size(E,1)];
  j = [E(:,1)' E(:,2)'];
  v = [ones(size(E,1),1)', -ones(size(E,1),1)'];
  % for repeated indices (i,j) sparse automatically sums up elements, as we want
  G = sparse(i,j,v,size(E,1),size(V,1));
end
