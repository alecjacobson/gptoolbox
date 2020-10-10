function [VV,QQ,SS,J] = dual_subdivide(V,F)
  % [VV,QQ,SS,J] = dual_subdivide(V,F)
  %
  % This is really just one iteration of Catmull Clark subdivision without
  % moving vertices.
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh triangle indices into V
  % Outputs:
  %   VV  #VV by 3 list of mesh vertex positions
  %   QQ  #QQ by 4 list of mesh quad indices into VV, QQ(:,1) are vertices of
  %     the original mesh (V).
  %   SS  #VV by #V matrix so that VV = SS*V
  %   J  #QQ list of indices of birth parents
  %   
  
  nv = size(V,1);
  nf = size(F,1);
  allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
  [uE,~,EMAP] = unique(sort(allE,2),'rows');
  EMAP = reshape(EMAP,size(F));
  ne = size(uE,1);
  SS = [ ...
    speye(nv,nv); ...
    sparse(repmat(1:ne,2,1)',uE,0.5,ne,nv); ...
    sparse(repmat(1:nf,3,1)',F,1/3,nf,nv); ...
    ];
  VV = SS*V;
  vec = @(X) reshape(X,[],1);
  J = repmat(1:nf,1,3)'; 
  I = nv+ne+J;
  QQ = [ ...
    vec(F(:,[1 2 3])) nv+vec(EMAP(:,[3 1 2])) I nv+vec(EMAP(:,[2 3 1]))];

end
