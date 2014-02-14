function B = barycentric_coordinates(P,V1,V2,V3)
  % BARYCENTRIC_COORDINATES Computes barycentric coordinates of point p in
  % simplex (v1,v2,v3)
  %
  % Inputs:
  %   P  #P by dim list of query point locations
  %   V1  #P by dim list of triangle corner locations
  %   V2  #P by dim list of triangle corner locations
  %   V3  #P by dim list of triangle corner locations
  % Outputs:
  %   B  #P by dim+1 list of barycentric coordinates
  %

  % SHOULD BE USING VARGIN TO ACCEPT ARBITRARY DIMENSION INPUT

  assert(size(P,1) == size(V1,1));
  assert(size(P,1) == size(V2,1));
  assert(size(P,1) == size(V3,1));
  % Only 2D supported for now
  assert(size(P,2) == 2);
  n = size(P,1);
  A1 = doublearea([ P;V2;V3],[1:n;n+[1:n;n+(1:n)]]');
  A2 = doublearea([V1; P;V3],[1:n;n+[1:n;n+(1:n)]]');
  A3 = doublearea([V1;V2; P],[1:n;n+[1:n;n+(1:n)]]');
  A  = doublearea([V1;V2;V3],[1:n;n+[1:n;n+(1:n)]]');
  B = bsxfun(@rdivide,[A1 A2 A3],A);
end
