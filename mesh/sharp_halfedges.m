function SHE = sharp_halfedges(V, F, angle, varargin)
% SHARP_HALFEDGES Compute sharp halfedges for a mesh based on dihedral angles.
%
% SHE = sharp_halfedges(V, F, angle)
% SHE = sharp_halfedges(V, F, angle, 'CornerIndexing')
% SHE = sharp_halfedges(V, F, angle, 'CyclicIndexing')
%
% This function marks "sharp" edges in a triangular mesh based on a given dihedral
% angle threshold. An edge is considered sharp if its dihedral angle exceeds the
% specified angle (in radians).
%
% Input:
%   V  #V by 3 matrix of vertex positions.
%   F  #F by 3 matrix of triangle vertex indices.
%   angle  Scalar threshold angle (in radians). An edge with a larger dihedral
%          angle is marked as sharp.
%
% Optional argument (mode):
%   'CornerIndexing' (default) - Uses corner-based indexing convention: edge j
%     is defined as the edge opposite vertex j in F(i,:).
%   'CyclicIndexing' - Uses edge-based cyclic indexing convention: local edge
%     indices follow the order of edges (1→2, 2→3, 3→1) in each triangle.
%
% Output:
%   SHE  #F by 3 matrix, where SHE(i,j) = 1 if the j-th edge of face i is sharp
%        (dihedral angle > angle), and 0 otherwise. Indexing depends on the mode.
%
% Example:
%   SHE = sharp_halfedges(V, F, pi/2-0.0001);
%   [EF, EI, uE, EMAP] = edge_flaps(F);
%   % isSharpEdge - mask of undirected sharp edges in uE.
%   isSharpEdge = zeros(size(uE,1),1);
%   isSharpEdge(EMAP(SHE==1)) = 1;

  % Default mode
  mode = 'CornerIndexing';

  % Parse optional argument
  if ~isempty(varargin)
    if ischar(varargin{1}) || isstring(varargin{1})
      mode = varargin{1};
    else
      error('Optional argument must be a string: ''CornerIndexing'' or ''CyclicIndexing''.');
    end
  end

  [A,C] = adjacency_dihedral_angle_matrix(V,F);
  [AI,AJ,AV] = find(A);
  keep = abs(AV-pi)>(angle) & ~isnan(AV);
  A = sparse(AI(keep),AJ(keep),1,size(A,1),size(A,2));
  [CI,~,CV] = find(C.*A);
  idx = sub2ind(size(C), CI, CV);
  SHE = zeros(size(F));
  SHE(idx)=1;

  % Apply cyclic reordering if requested
  if strcmpi(mode, 'CyclicIndexing')
    SHE = [SHE(:,3), SHE(:,1), SHE(:,2)];
  elseif ~strcmpi(mode, 'CornerIndexing')
    error('Unknown indexing mode. Use ''CornerIndexing'' or ''CyclicIndexing''.');
  end
end