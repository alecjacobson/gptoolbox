function T = flip_faces(T0,F)
  % T = flip_faces(T0,F)
  %
  % Inputs:
  %   T0  #T0 by 4 list of tet indices
  %   F  #F by 3 list of triangle face indices 
  % Outputs:
  %   T  #T' by 4 list of tet indices after flipping faces
  %

  T = T0;

  %    ...i...
  %  ./  / \  \.
  % a   /  .\---b
  %  \ /../  \ /
  %   j-------k
  % 
  % ijk lies on aijk and bikj
  %          ↓
  % abjk abik abji

  % edges before
  % ai
  % aj
  % ak
  % ij
  % ik
  % jk
  % bi
  % bj
  % bk
  %
  % edges after
  % ab
  % ai
  % aj
  % ak
  % bi
  % bj
  % bk
  % ij
  % ik
  % jk

  % O(T)
  allF = [ ...
    T(:,2) T(:,4) T(:,3); ...
    T(:,1) T(:,3) T(:,4); ...
    T(:,1) T(:,4) T(:,2); ...
    T(:,1) T(:,2) T(:,3); ...
    ];
  sF = sort(allF,2);
  [uF,~,FMAP] = unique(sF,'rows');
  A1 = allF(:,[1 2 3]);
  A2 = allF(:,[2 3 1]);
  A3 = allF(:,[3 1 2]);
  rev = ...
    ~(ismember(A1,uF,'rows') | ismember(A2,uF,'rows') | ismember(A3,uF,'rows'));

  F2T = sparse(FMAP,repmat(1:size(T,1),1,4)',T,size(uF,1),size(T,1));

  [found,to_flip] = ismember(sort(F,2),uF,'rows');
  assert(all(found),'All faces must be found in elements');

  [~,K] = find(F2T(to_flip,:));
  % no tets involved in multiple flips
  assert(numel(K) == numel(unique(K)),'Flipping faces involved in same tets not implemented');

  O = cumoccurrence(FMAP);
  opps = zeros(size(uF,1),2);
  opps(sub2ind(size(opps),FMAP,O)) = T(:);
  I = zeros(size(uF,1),2);
  I(sub2ind(size(I),FMAP,O)) = repmat(1:size(T,1),1,4)';
  J = zeros(size(uF,1),2);
  J(sub2ind(size(opps),FMAP,O)) = reshape(repmat(1:4,size(T,1),1),[],1);
  R = false(size(uF,1),2);
  R(sub2ind(size(opps),FMAP,O)) = rev;
  % <= because boundary edges only occur once (and might be correct
  assert(all(sum(R,2)<=1),'Must have consistent orientation');

  opps = opps(to_flip,:);
  I = I(to_flip,:);
  J = J(to_flip,:);
  R = R(to_flip,:);
  assert(all(I>0,'all'));
  IJKAB = [uF(to_flip,:) T(sub2ind(size(T),I,J))];
  IJKAB(R(:,1),[4 5]) = IJKAB(R(:,1),[5 4]);
  % New split tets
  ST = [IJKAB(:,[4 5 2 1]); ...
        IJKAB(:,[4 5 3 2]); ...
        IJKAB(:,[4 5 1 3])];
  T(I,:) = [];
  T = [T;ST];



end
