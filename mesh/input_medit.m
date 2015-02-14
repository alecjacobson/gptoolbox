function varargout = input_medit(V,F,C)
  % INPUT_MEDIT  Display input mesh using medit, highlighting boundary
  % edges and non-manifold edges
  %
  % C = input_medit(V,F,C)
  %
  % Inputs:
  %   V  #V by 3 list of mesh positions
  %   F  #F by 3 list of triangle indices
  %   Optional:
  %     C  #C connected component IDs
  %
  
  % Gather boundary and non-manifold edges
  NME = nonmanifold_edges(F);
  BE = outline(F);
  BE = setdiff(BE,NME,'rows');
  E = [NME 1+0*NME(:,1); ...
        BE 2+0*BE(:,1);];
  if ~exist('C') || isempty(C)
    C = randcycle(connected_components(F));
    C = C(:);
    % hack so that not all equal (medit doesn't like that)
    if all(C==C(1))
      C(ceil(rand*size(C,1))) = C(1)+1;
    end
  end
  if size(C,1)<size(V,1)
      C(size(V,1)) = 0;
      C(C==0) = max(C(:))+1;
  end
  % Gather self-intersecting faces
  ND = doublearea(V,F)>0;
  [~,~,IF] = selfintersect(V,F(ND,:),'DetectOnly',true);
  bad = false(size(F,1),1);
  bad(ND) = full(0<sparse(IF(:),1,1,sum(ND),1));
  F = [F bad];
  medit(V,[],F,'Data',C,'Edges',E,'Wait',false);
  if nargout>0
    varargout{1} = C;
  end
end
