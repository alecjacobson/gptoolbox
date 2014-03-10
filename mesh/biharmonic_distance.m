function [D] = biharmonic_distance(V,F,i,dim,p)
  % BIHARMONIC_DISTANCE Takes a mesh (V,F) and returns a distance field D from
  % all points in V to the ith vertex, according to the biharmonic embedding
  %
  % [D] = biharmonic_distance(V,F,i,dim)
  % [D] = biharmonic_distance(V,F,i,dim,p)
  %
  % Input:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices
  %   i  index of vertex from which to calculate distances
  %   dim  requested dimension of the embedding
  %   Optional:
  %     p  exponent above eigen values
  %       0.5  "semi-harmonic" embedding
  %       1  commute time embedding, "harmonic"
  %       2  biharmonic {default}
  %       3  "triharmonic" embedding
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

  % if power is not specfied use 2
  if(~exist('p','var'))
    p = 2;
  end

  B = biharmonic_embedding(V,F,dim,p);
  %B = biharmonic_embedding_yaron(V,F);
  D = sqrt(sum((repmat(B(i,:),size(B,1),1)-B(:,:)).^2,2));

  %tsurf(F,[V(:,1) V(:,2) D]);
  
end
