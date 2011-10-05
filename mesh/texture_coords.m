function [TC] = texture_coords(V,flip_y)
  % TEXTURE_COORDS normalize mesh x and y coordinates for use as texture
  % coordinates
  %
  % [TC] = texture_coords(V,flip_y)
  % 
  % Inputs:
  %   V  #V by 2 list of mesh vertex positions
  %   flip_y  boolean whether to flip texture coordinates, {true}
  % Outputs
  %   TC  #V by 2 list of mesh texture coordinates, between 0 and 1
  %
  % normalize texture coordinates to be between zero and 1
  %

  TC = [ ...
    (V(:,1)-min(V(:,1)))/(max(V(:,1))-min(V(:,1))) ...
    (V(:,2)-min(V(:,2)))/(max(V(:,2))-min(V(:,2)))];

  if ~exist('flip_y','var')
    flip_y = true;
  end

  if flip_y
    TC(:,2) = 1-TC(:,2);
  end
end
