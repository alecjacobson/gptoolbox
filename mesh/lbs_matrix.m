function M = lbs_matrix(V,W)
  % LBS_MATRIX  construct a matrix that when multiplied against a column of
  % affine transformation entries computes new coordinates of the vertices
  %
  % M = lbs_matrix(V,W)
  %
  % Input:
  %   V  #V by dim list of vertex rest positions
  %   W  #W by #handles list of correspondence weights
  % Output:
  %   M  #V * dim by #handles * dim * (dim+1) matrix such that
  %     new_V(:) = LBS(V,W,A) = reshape(M * A,size(V)), where A is a column
  %     vectors formed by the entries in each handle's dim by dim+1 
  %     transformation matrix. Specifcally, A =
  %       reshape(permute(Astack,[3 1 2]),n*dim*(dim+1),1)
  %     or A = [Lxx;Lyx;Lxy;Lyy;tx;ty], and likewise for other dim
  %     if Astack(:,:,i) is the dim by (dim+1) transformation at handle i
  %
  % Example:
  %  MLBS = lbs_matrix(V,W);
  %  n = size(V,1);
  %  dim = size(V,2);
  %  m = size(W,2);
  %  % stack of identity transformations
  %  Astack = repmat([eye(dim,dim) zeros(dim,1)],[1 1 m]);
  %  % collect transformations into column
  %  A = reshape(permute(Astack,[3 1 2]),m*dim*(dim+1),1);
  %  % apply transformations
  %  new_V = MLBS*A;
  %  new_V = reshape(new_V,[n dim]);
  %  
  %  % Alternative:
  %  % Q is #W by 4 list of quats, T is #W by 3 list of translations
  %  % A is #W*4 by 3 stack of transposed affine transformations
  %  A = reshape(cat(2,permute(quat2mat(Q),[2 1 3]),permute(T,[2 3 1])),3,[])';
  %  % M is #V by #W*4 skinning matrix
  %  M = reshape(bsxfun(@times,[V ones(n,1)],permute(W,[1 3 2])),n,[]);
  %  U = M*A;
  %  

  % number of mesh (domain) vertices
  n = size(V,1);
  assert(n == size(W,1));
  % dimension of mesh
  dim = size(V,2);
  % number of handles
  m = size(W,2);

  % M = zeros(V*dim,m*dim*(dim+1));

  % repeat vertex positions so that VV(:,:,i) gives #V by #handles matrix where
  % each column is ith coordinates of vertices
  VV = permute(repmat(V,[1 1 m]),[1 3 2]);
  % multiply each column in VV by respective weights for that handle
  VVW = VV.*repmat(W,[1 1 dim]);
  % matrix of zeros
  Z = zeros(n,m);
  switch dim
  case 2
    M = [ ...
      VVW(:,:,1) Z          VVW(:,:,2) Z          W Z; ...
      Z          VVW(:,:,1) Z          VVW(:,:,2) Z W];
  case 3
    M = [ ...
      VVW(:,:,1) Z          Z          VVW(:,:,2) Z          Z          VVW(:,:,3) Z          Z          W Z Z; ...
      Z          VVW(:,:,1) Z          Z          VVW(:,:,2) Z          Z          VVW(:,:,3) Z          Z W Z; ...
      Z          Z          VVW(:,:,1) Z          Z          VVW(:,:,2) Z          Z          VVW(:,:,3) Z Z W];
  otherwise
    error('Only dim=2 or dim=3 supported');
  end


end
