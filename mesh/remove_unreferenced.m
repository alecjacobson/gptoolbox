function [RV,IM] = remove_unreferenced(V,F)
  % REMOVE_UNREFERENCED Removes any rows in V that are not referenced in R
  %
  % [RV,IM] = remove_unreferenced(V,F)
  %
  % Inputs:
  %   V  #V by dim list of "vertex positions"
  %   F  #F by anything list of indices into V (will be treated as F(:))
  % Outputs:
  %  RV  #V by dim vertex positions, order such that if the jth vertex is
  %    some face in F, and the kth vertex is not then j comes before k
  %  IM  #V by 1 list of indices such that: RF = IM(F) and RT = IM(T)
  %    and V(find(IM<=size(SV,1)),:) = RV
  % 
  % Examples:
  %   % Tet mesh in (V,T,F)
  %   [RV,I] = remove_unreferenced(V,[T(:);F(:)]);
  %   T = I(T);
  %   F = I(F);
  %   ... % do some computation on RV
  %   % replace back into V
  %   V(find(IM<=size(SV,1)),:) = V
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
  % reorder vertices
  RV(IM,:) = V;
  % Remove unreferenced
  RV = RV(1:max(IM(F(:))),:);
end
