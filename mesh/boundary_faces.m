function [F,J] = boundary_faces(T)
  % BOUNDARY_FACES Determine boundary faces of tetrahedra stored in T
  %
  % F = boundary_faces(T)
  %
  % Input:
  %   T  tetrahedron index list, m by 4, where m is the number of tetrahedra
  % Output:
  %   F  list of boundary faces, n by 3, where n is the number of boundary faces
  %   J  list of indices into T, n by 1
  %

  % get all faces
  allF = [ ...
    T(:,1) T(:,2) T(:,3); ...
    T(:,1) T(:,3) T(:,4); ...
    T(:,1) T(:,4) T(:,2); ...
    T(:,2) T(:,4) T(:,3)];
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
  J = mod(sorted_exteriorJ(J(I))-1,size(T,1))+1;
  F = allF(I,:);
end
