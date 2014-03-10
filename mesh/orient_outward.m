function [FF,I] = orient_outward(V,F,C)
  % ORIENT_OUTWARD Use a simple heuristic to maintain or flip each
  % manifold/orientable patch of a mesh. Assumes that independent patches have
  % already been oriented consistently up to sign.
  % 
  % [FF,I] = orient_outward(V,F,C)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  %   C  #F list of component ids
  % Outputs:
  %   FF  #F by 3 list of new, potentially flipped triangle indices
  %   I  #C list of bools whether patch was flipped
  %
  % See also: bfs_orient, manifold_patches
  %
  [FF,C] = bfs_orient(F);
  I = false(max(C),1);
  for c = 1:max(C)
    N = normalizerow(normals(V,FF(C==c,:)));
    A = doublearea(V,FF(C==c,:));
    BC = barycenter(V,FF(C==c,:));
    BCmean = A'*BC/sum(A);
    BC = bsxfun(@minus,BC,BCmean);
    ndot = sum(bsxfun(@times,A,sum(N.*BC,2)));
    if ndot<0
      FF(C==c,:) = fliplr(FF(C==c,:));
      I(c) = true;
    else
      I(c) = false;
    end
  end
end
