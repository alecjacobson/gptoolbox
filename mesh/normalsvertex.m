function [ N ] = normalsvertex( V, F )
  % NORMALSVERTEX 
  %
  % N = normalsvertex(V,F)
  %
  % Compute normals per vertex with angle weights
  %
  % Inputs:
  %  V  #V x 3 matrix of vertex coordinates
  %  F  #F x 3  matrix of indices of triangle corners
  % Output:
  %  N  #F x 3 list of vertex normals
  %
  
  i1 = F(:,1);
  i2 = F(:,2);
  i3 = F(:,3);
  
  Nf = normals(V,F);
  Nf = normalizerow(Nf);
  
  N = zeros(size(V));
  
  A = internalangles(V,F);
  
  for i=1:size(F,1) % for each face
      Fi = F(i,:);
      
      N(Fi(1),:) = N(Fi(1),:) + (Nf(i,:) * A(i));
      N(Fi(2),:) = N(Fi(2),:) + (Nf(i,:) * A(i));
      N(Fi(3),:) = N(Fi(3),:) + (Nf(i,:) * A(i));
      
  end
  
end

