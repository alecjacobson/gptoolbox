function [offsets,WI,parents,I] = bone_forest(C,E,P,BE)
  % [offsets,WI,parents] = bone_forest(C,E,P,BE)
  %
  % See also: writeBF.m

  if nargin == 2 || isempty(BE)
    BE = E;
    P = [];
  end
  [BE,I] = sortrows(BE);
  % point properties
  P_parents = zeros(numel(P),1);
  P_offsets = C(P,:);
  P_WI = (1:numel(P))';
  % bone properties
  % bones form tree(s)
  assert(numel(unique(BE(:,2))) == numel(BE(:,2)));
  % 0 if not dest otherwise index in list
  index_map = zeros(max(BE(:)),1);
  index_map(BE(:,2)) = 1:size(BE,1);
  roots = setdiff(BE(:,1),BE(:,2));
  % add roots to index map
  index_map(roots) = size(BE,1)+(1:numel(roots));
  % parents after adding roots
  BE_parents = index_map(BE(:,1));
  % roots have no mother or father
  BE_parents = [BE_parents; zeros(numel(roots),1)];
  BE_WI = [1:size(BE,1) zeros(1,numel(roots))]';
  BE_offsets = bone_offsets(C([BE(:,2); roots],:),BE_parents);
  % we'll first put "point" handles
  % then bone's (recall that in BF each joint needs an entry, including roots)
  parents = [P_parents;numel(P).*(BE_WI~=0)+BE_parents];
  WI = [P_WI;numel(P).*(BE_WI~=0)+BE_WI];
  offsets = [P_offsets;BE_offsets];
end
