function [C,V,D] = commute_time_embedding(V,F,dim);
  % [C,V,D] = commute_time_embedding(V,F);
  %
  % Takes a mesh (V,F) and returns an embedding using the spectrum of the
  % harmonic operator. Then the commute time distance between two points i and
  % j can be computed as the euclidean distance between C(i,:) and C(j,:),
  % namely: 
  %   dist_ij = sqrt(sum((C(i,:)-C(j,:)).^2,2));
  % We also once referred to this as "harmonic embedding" because "Biharmonic
  % Distance" uses an embedding where the eigenvalues in the denominator are
  % squared, hence that was the "biharmonic embedding", thus this could be the
  % "harmonic embedding". Not to be confused with the embedding used in
  % diffusion distance.
  % 
  % Input:
  %   V  vertex list
  %   F  face list
  %   dim requested dimension of the embedding
  % Output:
  %   C  commute time embedding
  %   V  eigenvectors used in embedding
  %   D  eigenvalues used in embedding
  % 

  % if dimension is not specfied use 4
  if(~exist('dim','var'))
    dim = 4;
  end

  % get cotangent matrix
  L = cotmatrix(V,F);
  % get mass matrix
  M = massmatrix(V,F,'voronoi');
  % get dim+1 smallest magnitude eigenvalues and corresponding vectors
  %[V,D] = eigs(L,M,dim+1,'sm');
  %V = V(:, 2:end);
  %D = D(2:end, 2:end);

  % This also works, because of the sign change in the eigenvalues matlab
  % reverses the output order so 0.0 is the last eigenvalue
  [V,D] = eigs(-2*L,M./sum(M(:)),dim+1,'sm');
  V = V(:, 1:end-1);
  D = D(1:end-1, 1:end-1);

  %  divide each eigenvector by corresponding eigenvalue
  C = V * (inv(D));
end
