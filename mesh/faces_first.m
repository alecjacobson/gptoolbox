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
  %   V  #vertices by 3 vertex positions
  %   T  #tetrahedra by 4 list of tetrahedra indices
  %   F  #faces by 3 list of face indices
  % Output: 
  %   RV  #vertices by 3 vertex positions, order such that if the jth vertex is
  %     some face in F, and the kth vertex is not then j comes before k
  %   RT  #tetrahedra by 4 list of tetrahedra indices, reindexed to use RV
  %   RF  #faces by 3 list of face indices, reindexed to use RV
  %   IM  #vertices by 1 list of indices such that: RF = IM(F) and RT = IM(T)
  %     and RV(IM,:) = V
  %

  [~,IM] = remove_unreferenced(V,F(:));

  % reorder vertices
  RV(IM,:) = V;
  % reindex faces: reshape is necessary for #F = 1
  RF = reshape(IM(F),size(F));
  % reindex tetrahedra
  RT = reshape(IM(T),size(T));
end
