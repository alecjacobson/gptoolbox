function Z = iminpaint(im,M)
% Z = iminpaint(im,M)
% 
% Inputs:
%   im  w by h by c image
%   M  w by h mask
% Outputs:
%   Z  w by h by c image with mask in-painted with Laplacian
%
  if islogical(M)
    assert(numel(M) == size(im,1)*size(im,2));
    if any(size(M)==1)
      M = reshape(M,size(im,1),size(im,2));
    end
  end
  L = fd_laplacian(size(im(:,:,1)));
  b = find(M);
  X = reshape(im,[],size(im,3));
  Z = min_quad_with_fixed(-L,[],b,X(b,:),[],[],struct('force_Aeq_li',true));
  Z = reshape(Z,size(im));
end
