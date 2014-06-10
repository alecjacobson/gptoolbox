function [ B ] = barycenter(V, F)
  % BARYCENTER Compute the barycenter of every triangle
  %
  % B = barycenter(V,F)
  %
  % Inputs:
  %   V #V x dim matrix of vertex coordinates
  %   F #F x simplex_size  matrix of indices of triangle corners
  % Output:
  %   B a #F x dim matrix of 3d vertices
  % 
  % See also: quadrature_points
  %

  B = zeros(size(F,1),size(V,2));
  for ii = 1:size(F,2)
    B = B + 1/size(F,2) * V(F(:,ii),:);
  end

end

