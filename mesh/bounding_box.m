function [BB,BF] = bounding_box(V)
  % BOUNDING_BOX  Compute the bounding box of a set of points
  % 
  % [BB,BF] = bounding_box(V)
  %
  % Inputs:
  %   V  #V by dim list of points
  % Outputs:
  %   BB 2^dim list of boundary box vertices
  %   BF #BF by #dim list of facets
  %

  dim = size(V,2);
  minV = min(V);
  maxV = max(V);
  C = combn([0 1],dim) == 1;
  BB = bsxfun(@times,C,minV) + bsxfun(@times,~C,maxV);
  
  % lazy faces
  D = DelaunayTri(BB);
  switch size(D.Triangulation,2)
  case 3
    BF = outline(D.Triangulation);
  case 4
    BF = boundary_faces(D.Triangulation);
  otherwise
    error('Not supported');
  end



end
