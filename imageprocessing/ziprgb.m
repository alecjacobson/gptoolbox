function rgb = ziprgb(R,G,B)
  % Takes R,G,B as nxm matrices and makes rgb a nxmx3 3d-matrix
  rgb(:,:,1) = R;rgb(:,:,2) = G;rgb(:,:,3) = B;
end
