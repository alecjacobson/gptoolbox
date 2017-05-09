function L = fd_laplacian(side)
  % FD_LAPLACIAN  build a finite difference laplacian for a regular grid.
  %
  % L = fd_laplacian([h,w])
  % L = fd_laplacian([h,w,t])
  % 
  % Inputs:
  %   dims  number of nodes along height (and width and depth)
  % Outputs:
  %   L   prod(side) by prod(side) Laplacian with negative diagonal.
  % 
  % See also: fd_grad

  function B = vec(A)
    B = A(:)';
  end

  switch numel(side)
  case 2
    h = side(1);
    w = side(2);
    I = bsxfun(@plus,h*(0:w-1),(1:h)');
    L = sparse( ...
      [vec(I(2:end,:)) vec(I(:,2:end))], ...
      [vec(I(1:end-1,:)) vec(I(:,1:end-1))], ...
      1,w*h,w*h);
  case 3
    I = sub2ind(side([2 1 3]),1:prod(side))';
    I = reshape(I,side);
    L = sparse( ...
      [vec(I(2:end  ,:,:)) vec(I(:,2:end  ,:)) vec(I(:,:,2:end  ))], ...
      [vec(I(1:end-1,:,:)) vec(I(:,1:end-1,:)) vec(I(:,:,1:end-1))], ...
      1,prod(side),prod(side));
  end
  L = L + L';
  L = L-diag(sum(L,2));
end
