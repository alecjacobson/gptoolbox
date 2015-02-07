function [F,J,K] = boundary_faces(T)
  % BOUNDARY_FACES Determine boundary faces of tetrahedra stored in T
  %
  % F = boundary_faces(T)
  %
  % Input:
  %   T  tetrahedron index list, m by 4, where m is the number of tetrahedra
  % Output:
  %   F  list of boundary faces, n by 3, where n is the number of boundary faces
  %   J  list of indices into T, n by 1
  %   K  list of indices revealing across from which vertex is this face
  %

  % get all faces
  allF = [ ...
    T(:,[4 2 3]); ...
    T(:,[3 1 4]); ...
    T(:,[2 4 1]); ...
    T(:,[1 3 2])];
  % sort rows so that faces are reorder in ascending order of indices
  sortedF = sort(allF,2);
  % determine uniqueness of faces
  [u,m,n] = unique(sortedF,'rows');
  % determine counts for each unique face
  counts = accumarray(n(:), 1);
  % extract faces that only occurred once
  sorted_exteriorF = u(counts == 1,:);
  sorted_exteriorJ = m(counts == 1,:);
  % find in original faces so that ordering of indices is correct
  [I,J] = ismember(sortedF,sorted_exteriorF,'rows');
  L = sorted_exteriorJ(J(I))-1;
  J = mod(L,size(T,1))+1;
  K = floor(L/size(T,1))+1;
  F = allF(I,:);
end
