function P = bone_parents(BE)
  % BONE_PARENTS Find parents in bone
  %
  % P = bone_parents(BE)
  %
  % Inputs:
  %   BE  #BE by 2 list of **directed** bone edge indices 
  % Outputs:
  %   P  #BE list of parent indices into BE (0 means root)
  %
  %

  %BE = sortrows(BE);
  % 0 if not dest otherwise index in list
  index_map = zeros(max(BE(:)),1);
  index_map(BE(:,2)) = 1:size(BE,1);
  roots = setdiff(BE(:,1),BE(:,2));
  % add roots to index map
  index_map(roots) = 0;%size(BE,1)+(1:numel(roots));
  % parents after adding roots
  P = index_map(BE(:,1));
end
