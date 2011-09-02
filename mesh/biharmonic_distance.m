function [D] = biharmonic_distance(V,F,i,dim)
  % [D] = biharmonic_distance(V,F,i,dim)
  %
  % Takes a mesh (V,F) and returns a distance field D from all points in V to
  % the ith vertex, according to the biharmonic embedding
  %
  % Input:
  %   V  vertex list
  %   F  face list
  %   i  index of vertex from which to calculate distances
  %   dim  requested dimension of the embedding
  % Output:
  %   D  biharmonic distance field 
  % 

  % if index not given then use n/2th point
  if(~exist('i','var'))
    i = ceil(size(V,1)/2 + sqrt(size(V,1))/2);
  end
  
  % if dimension is not specfied use 4
  if(~exist('dim','var'))
    dim = 4;
  end

  B = biharmonic_embedding(V,F,dim);
  D = sqrt(sum((repmat(B(i,:),size(B,1),1)-B(:,:)).^2,2));

  tsurf(F,[V(:,1) V(:,2) D]);
  
end
