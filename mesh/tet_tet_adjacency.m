function [Tp,Ti,uF,FMAP] = tet_tet_adjacency(T)
  % [Tp,Ti,uF,FMAP] = tet_tet_adjacency(T);
  %
  % Inputs:
  %   T  #T by 4 list of tet indices into rows of some V
  % Outputs:
  %   Tp  #T by 4 list of tet indices of adjacent tets (0 if boundary)
  %   Ti  #T by 4 list of local face indices of adjacent tets (0 if boundary)
  %   uF  #F by 3 list of unique tet face indices into rows of some V
  %   FMAP  #T*4 by 1 list of indices into uF for each face of T
  %
  % See also: crouzeix_raviart_cotmatrix
  %
  F = [ ...
    T(:,2) T(:,4) T(:,3); ...
    T(:,1) T(:,3) T(:,4); ...
    T(:,1) T(:,4) T(:,2); ...
    T(:,1) T(:,2) T(:,3); ...
    ];
  % Map duplicate facets to first instance
  %[F,~,FMAP] = unique(sort(allF,2),'rows');
  [uF,~,FMAP] = unique(sort(F,2),'rows');
  K = reshape(repmat(1:4,size(T,1),1),4*size(T,1),1);
  T2F = sparse(FMAP,repmat(1:size(T,1),1,4)',K);
  [I1,J1] = max(T2F,[],2);
  T2F(sub2ind(size(T2F),(1:size(uF,1))',J1)) = 0;
  [I2,J2] = max(T2F,[],2);
  int = find(I2>0);
  J1 = J1(int);
  J2 = J2(int);
  I1 = I1(int);
  I2 = I2(int);
  Tp = zeros(size(T));
  Ti = zeros(size(T));
  Tp(sub2ind(size(Tp),J1,I1)) = J2;
  Tp(sub2ind(size(Tp),J2,I2)) = J1;
  Ti(sub2ind(size(Ti),J1,I1)) = I2;
  Ti(sub2ind(size(Ti),J2,I2)) = I1;
end
