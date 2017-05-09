function [L,NUM] = gp_bwlabel(BW,N)
  % GP_BWLABEL Label connected components in the true regions of a logical
  % BWage. Assuming [0 1 0;1 1 1;0 1 0] connectivity.
  %
  % L = gp_bwlabel(BW)
  %
  % Inputs:
  %   BW  h by w logical BWage
  % Outputs:
  %   L  h by w BWage of ids of regions (0 for "unlabeled")
  %
  % See also: gp_bwboundaries
  %
  assert(islogical(BW));

  if nargin<2
    N = 4;
  end

  % Adjacency matrix
  A = fd_laplacian(size(BW))~=0;
  if N == 8
    A = (A*1*A)>1;
  end
  A(~BW,:) = 0;
  A(:,~BW) = 0;
  [~,L] = conncomp(A);
  assert(min(L(:)) > 0);
  L(~BW(:)) = 0;
  [~,~,L] = unique(L);
  L = reshape(L-1,size(BW));
  assert(all(L(~BW(:))==0));
  NUM = max(L(:));
end
