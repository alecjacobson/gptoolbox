function A = vector_area_matrix(F)
  % Constructs the symmetric area matrix A, s.t.  V(:)' * A * V(:) is the
  % **vector area** of the mesh (V,F).
  %
  % A = vector_area_matrix(F)
  %
  % Inputs:
  %   F  #F by 3 list of mesh faces (must be triangles)
  % Outputs:
  %   A  #Vx2 by #Vx2 area matrix
  %
  % number of vertices
  n = max(F(:));
  O = outline(F);
  A = sparse( ...
    [O;O(:,[2 1])+n],[O(:,[2 1])+n;O],repmat([1 -1]/4,size(O,1)*2,1),2*n,2*n);
end
