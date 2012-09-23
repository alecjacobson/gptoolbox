function [RV,RT,RF,IM] = faces_first(V,T,F)
  % FACES_FIRST Reorder vertices so that vertices in face list come before
  % vertices that don't appear in the face list. This is especially useful if
  % the face list contains only surface faces and you want surface vertices
  % listed before internal vertices
  %
  % [RV,RT,RF,IM] = faces_first(V,T,F);
  %
  %
  % Input:
  %  V  # vertices by 3 vertex positions
  %  T  # tetrahedra by 4 list of tetrahedra indices
  %  F  # faces by 3 list of face indices
  % Output: 
  %  RV  # vertices by 3 vertex positions, order such that if the jth vertex is
  %    some face in F, and the kth vertex is not then j comes before k
  %  RT  # tetrahedra by 4 list of tetrahedra indices, reindexed to use RV
  %  RF  # faces by 3 list of face indices, reindexed to use RV
  %  IM  # vertices by 1 list of indices such that: RF = IM(F) and RT = IM(T)
  %    and RV(IM,:) = V
  %

  % get list of unique vertex indices that occur in faces
  U = unique(F(:));
  % get list of vertices that do not occur in faces
  NU = (1:size(V,1))';
  NU = NU(~ismember(NU,U));
  assert((size(U,1) + size(NU,1)) == size(V,1));
  % allocate space for an indexmap so that IM[i] gives new index of vertex i
  IM = zeros(size(V,1),1);
  % reindex vertices that occur in faces to be first
  IM(U) = 1:size(U,1);
  % reindex vertices that do not occur in faces to come after those that do
  IM(NU) = size(U,1) + (1:size(NU,1));
  % reindex faces
  RF = IM(F);
  % reindex tetrahedra
  RT = IM(T);
  % reorder vertices
  RV(IM,:) = V;
end
