function [VV,QQ,SS,J] = catmull_clark(V,Q,iter)
  % CATMULL_CLARK Perform iter iterations of catmull-clark subdivision on a
  % regular, manifold quad mesh (V,Q)
  %
  % Inputs:
  %   V  #V by dim list of mesh vertex positions
  %   Q  #Q by 4 list of quad indices into V
  % Outputs:
  %   VV  #VV by dim list of output mesh vertex positions
  %   QQ  #QQ by 4 list of quad mesh indices into VV
  %   SS  #VV by #V linear subdivision operator so that VV = SS*V
  %   J  #QQ list of indices of birth parents into Q
  % 

  if (~exist('iter','var'))
      iter = 1;
  end
  VV = V;
  SS = speye(size(V,1));
  QQ = Q;
  J = (1:size(Q,1))';

  for i=1:iter
    % assume all are quads OR tris
    FF = QQ;
      
    % number of original vertices
    n = size(VV,1);
    % number of original faces
    nf = size(FF,1);

    % extract pure quads
    pure = false(nf,1);
    if size(FF,2)==4
      pure = FF(:,4)~=nan & FF(:,4)>0 & FF(:,4)~=FF(:,3);
      QQ = FF(pure,:);
    else
      QQ = zeros(0,4);
    end
    TT = FF(~pure,1:3);
    JQ = J(pure);
    JT = J(~pure);
    nq = size(QQ,1);
    nt = size(TT,1);
    assert(nf == (nq+nt));

    % get "quad edges" aka half-edges
    QE = [QQ(:,1:2);QQ(:,2:3);QQ(:,3:4);QQ(:,[4 1])];
    % get "triangle edges"
    TE = [TT(:,1:2);TT(:,2:3);TT(:,[3 1])];
    % all edges
    FE = [QE;TE];
    % get unique edge list, QE2E tells original face edges (with duplicates,
    % EE) where to find matching unique edge (E)
    [E,I,FE2E] = unique(sort(FE,2),'rows');
    FEF = [repmat(1:nq,1,4) nq+repmat(1:nt,1,3)]';
    % number of original edges
    ne = size(E,1);
    E2F = zeros(ne,2);
    flip = FE(:,1)~=E(FE2E,1);
    E2F(sub2ind(size(E2F),FE2E,flip+1)) = FEF;

    % face barycenter
    SF = [ ...
      sparse(repmat(1:nq,1,4),QQ(:)',1/4*ones(1,nq*4),nq,n); ...
      sparse(repmat(1:nt,1,3),TT(:)',1/3*ones(1,nt*3),nt,n)];
    % Average of original edge points and edge-flap face barycenters
    SE = sparse(repmat(1:ne,1,4),[E(:);n+E2F(:)]',1/4*ones(1,ne*4),ne,n+nf) ...
      *[speye(n);SF];

    V2F = [ ...
      sparse(QQ,repmat(1:nq,4,1)',1,n,nq) ...
      sparse(TT,repmat(1:nt,3,1)',1,n,nt)];
    % number of faces incident on each vertex
    val = sum(V2F,2);
    % normalize to take average
    V2F = bsxfun(@rdivide,V2F,sum(V2F,2));
    % compute midpoints of original vertices
    ME = full(sparse(E(:,[1 2 1 2]),E(:,[1 2 2 1]),1,n,n));
    ME = bsxfun(@rdivide,ME,sum(ME,2));
    % Weighting 
    W = bsxfun(@rdivide,[(val-3) ones(n,1) 2*ones(n,1)],val);
    SV = ...
      sparse(repmat(1:n,3,1)',reshape((1:3*n)',n,3), W,n,3*n)* ...
      [speye(n);V2F*SF;ME];

    S = [SV;SF;SE];
    VV = S*VV;

    %clf;
    %hold on;
    %tsurf(Q,V,'FaceAlpha',0.5);
    %scatter3(VV(:,1),VV(:,2),VV(:,3),'LineWidth',3);
    %hold off;
    %error

    SS = S*SS;
    QE2E = reshape(FE2E(1:nq*4),size(QQ));
    TE2E = reshape(FE2E(nq*4+(1:nt*3)),size(TT));
    QQ = [ ...
      QQ(:,1) n+nf+QE2E(:,1) n+(1:nq)' n+nf+QE2E(:,4); ...
      QQ(:,2) n+nf+QE2E(:,2) n+(1:nq)' n+nf+QE2E(:,1); ...
      QQ(:,3) n+nf+QE2E(:,3) n+(1:nq)' n+nf+QE2E(:,2); ...
      QQ(:,4) n+nf+QE2E(:,4) n+(1:nq)' n+nf+QE2E(:,3); ...
      ];
    TT = [ ...
      TT(:,1) n+nf+TE2E(:,1) n+nq+(1:nt)' n+nf+TE2E(:,3); ...
      TT(:,2) n+nf+TE2E(:,2) n+nq+(1:nt)' n+nf+TE2E(:,1); ...
      TT(:,3) n+nf+TE2E(:,3) n+nq+(1:nt)' n+nf+TE2E(:,2)];
    QQ = [QQ;TT];
    J = [JQ;JQ;JQ;JQ;JT;JT;JT];

  end
      
    
end
