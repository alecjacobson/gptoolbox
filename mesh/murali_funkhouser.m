function [S,DV,DEle,DF,DN] = murali_funkhouser(V,F,varargin)
  % Compute "solidness" according to "Consistent solid and boundary
  % representations from arbitrary polygonal data" [Murali et Funkhouser 1997]
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by dim list of facets indexing V
  %   Optional:
  %     'AlreadyClean' followed by true of false. Only applies to dim==3.
  % Outputs:
  %   S  #DEle list of scalar function values
  %   DV  #DV by dim list of space decomposition vertices
  %   DEle  #DEle by dim+1 list of space decomposition elements indexing DV
  %   DF  #DF by dim list of facets (whose union forms input F)
  %   DN  #DEle by dim+1 list of space decomposition neighbors indexing DEle
  % 
  % Example:
  %   [S,DEle,DV,DF] = murali_funkhouser(V,F);
  %   % visualize output
  %   medit(DV,DEle,DF,'Data',S);
  %   % "Solid" tets
  %   CEle = DEle(S>0,:);
  %   [CV,I] = remove_unreferenced(DV,CEle);
  %   CEle = I(CEle);
  %   % Accordingly oriented boundary faces
  %   CF = boundary_faces(CEle);
  %   % Visualize solution
  %   medit(CV,CEle,CF);
  %   
  %

  already_clean = false;
  v = 1;
  while v <= numel(varargin)
    switch varargin{v}
    case 'AlreadyClean'
      assert((v+1)<=numel(varargin));
      v = v+1;
      already_clean = v;
    otherwise
      error(sprintf('Unsupported input parameter "%s"',varargin{v}));
    end
    v = v+1;
  end

  % simplex size
  dim = size(F,2);
  % Facet areas
  switch dim
  case 3
    if already_clean
      SV = V;
      SF = F;
    else
      % print some statistics
      statistics(V,F)
      % Clean up mesh for TetGen
      [SV,SF,SVJ] = clean(V,F, ...
        'MinDist',1e-6,'MinAngle',-1,'MinArea',-1,'SelfIntersections','mesh');
    end
    % Use TetGen to mesh CDT
    [DV,DEle,DF,DN] = cdt(SV,SF,'TetgenFlags','-Y','UseBoundingBox',true,'BoundingBoxUpsample',10,'BoundingBoxPush',3);
    allF = [ ...
      DEle(:,[2 3 4]); ...
      DEle(:,[3 4 1]); ...
      DEle(:,[4 1 2]); ...
      DEle(:,[1 2 3])];
    A = reshape(doublearea(DV,allF),size(DEle));
  case 2
    %[DV,DEle,DF,DN] = cdt(V,F,'TriangleFlags','-Y','UseBoundingBox',true);
    [DV,DEle,DF,DN] = cdt(V,F,'TriangleFlags','-q31','UseBoundingBox',true, ...
      'BoundingBoxPush',1.1);
    [DV,DEle,DF,DN] = cdt(V,F, ...
      'TriangleFlags',sprintf('-q32 -a%f',mean(doublearea(DV,DEle))/8), ...
      'UseBoundingBox',true, ...
      'BoundingBoxPush',1.1);
    allF = [ ...
      DEle(:,[2 3]); ...
      DEle(:,[3 1]); ...
      DEle(:,[1 2]);];
    A = reshape(sqrt(sum((DV(allF(:,1),:)-DV(allF(:,2),:)).^2,2)),size(DEle));
  end
  % Union of DF == SF
  [I] = in_elements(DF,DEle);
  assert(all(I));
  % Which faces of allF are in SF
  inSF = ismember(sort(allF,2),sort(DF,2),'rows');
  % si = (∑j (tij - oij) sj) / Ai
  % 0 = (∑j (tij - oij) sj) / Ai - si
  % 0 = (∑j (tij - oij) sj - aij si) / Ai
  % : tij + oij = aij
  %   tij - oij = aij - 2*oij
  % 0 = (∑j (aij - 2*oij) sj - aij si) / Ai
  % 0 = (∑j (aij - 2*oij) sj - (aij - 2*oij + 2*oij) si) / Ai
  % 0 = (∑j (aij - 2*oij) sj - (aij - 2*oij) si - (2*oij) si) / Ai
  % 0 = (∑j (aij - 2*oij) (sj - si) - 2*oij si) / Ai
  % 0 = (∑j (aij - 2*oij) (sj - si))/ Ai  - si ∑ 2*oij / Ai
  % 0 = (∑j (tij - oij) (sj - si))/ Ai  - si ∑ 2*oij / Ai
  %
  % (L - D) * s = 0
  % 
  % Where L is area-weighted laplacian and D is a area-weighted diagonal
  %
  % 0 = 1/Ai * ( (∑j (tij - oij) (sj - si))  - si ∑ 2*oij)
  %
  % IA * (L - D) * s = 0
  % 
  % Where L is area-weighted laplacian and D is a area-weighted diagonal, IA is
  % a diagonal matrix with 1/Ai.
  %
  % So it's like a Laplacian minus a positive diagonal.
  %warning('uniform');
  %A(:) = 1;
  T = reshape(A,size(DEle,1)*(dim+1),1);
  O = reshape(A,size(DEle,1)*(dim+1),1);
  T(inSF) = 0;
  O(~inSF) = 0;
  Ai = A;
  Ai(DN==-1) = 0;
  Ai = sum(Ai,2);
  % Ai si + (∑j (oij - tij) sj) = 0
  % Mii = Ai
  % Mij = oij - tij
  % off-diagonal then diagonal
  MI = repmat(1:size(DEle,1),1,(dim+1)+1)';
  MJ = [DN(:);(1:size(DEle,1))'];
  MV = [O(:) - T(:);Ai(:)];
  %MV = [-A(:);Ai(:)];
  M = sparse(MI(MJ>0),MJ(MJ>0),MV(MJ>0), size(DEle,1),size(DEle,1));
  assert(max(max(abs(M-M')))<1e-10);

  % Fix boundary elements whose boundary facet is not original to outside
  [~,BC] = on_boundary(DEle);
  b = find(any(BC>0,2) & ~any(BC>0 & reshape(inSF,size(DEle)),2));
  in = setdiff(1:size(DEle,1),b);
  % Solve 
  [L,p,P] = chol(M(in,in),'lower');
  S = zeros(size(DEle,1),1);
  S(b) = -1;
  S(in) = P * (L' \ (L \ ( P' * (-M(in,b) * S(b)))));

  %switch dim
  %case 3
  %  medit(DV,DEle,boundary_faces(DEle),'Data',S)
  %case 2
  % set(tsurf(DEle,DV),'CData',S);
  % %colormap(jet(ceil(max(S)-min(S))));
  % colormap(jet(256));
  % axis equal
  % colorbar
  %end

end
