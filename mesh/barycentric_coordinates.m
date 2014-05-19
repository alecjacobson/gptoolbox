function B = barycentric_coordinates(P,varargin)
  % BARYCENTRIC_COORDINATES Computes barycentric coordinates of point p in
  % simplex (v1,v2,v3)
  % 
  % B = barycentric_coordinates(P,V1,V2,V3, ...)
  %
  % Inputs:
  %   P  #P by dim list of query point locations
  %   V1  #P by dim list of simplex corner locations
  %   V2  #P by dim list of simplex corner locations
  %   ...
  % Outputs:
  %   B  #P by dim+1 list of barycentric coordinates
  %

  % SHOULD BE USING VARGIN TO ACCEPT ARBITRARY DIMENSION INPUT
  for varg = varargin
    assert(size(P,1) == size(varg{1},1), 'All inputs should be same length');
    assert(size(P,2) == size(varg{1},2), 'All inputs should be same dimension');
  end

  n = size(P,1);
  switch numel(varargin)
  case 4
    V1 = varargin{1};
    V2 = varargin{2};
    V3 = varargin{3};
    V4 = varargin{4};
    T = bsxfun(@plus,(1:n)',(0:3)*n);
    A1 = volume([V2;V4;V3;P],T);
    A2 = volume([V1;V3;V4;P],T);
    A3 = volume([V1;V4;V2;P],T);
    A4 = volume([V1;V2;V3;P],T);
    A  = volume([V1;V2;V3;V4],T);
    if size(P,2)>3 && max(abs(sum([A1 A2 A3 A4],2)-A))>1e-14
      warning('Possibly negative coordinates. Not supported in dim~=3');
    end
    B = bsxfun(@rdivide,[A1 A2 A3 A4],A);
  case 3
    V1 = varargin{1};
    V2 = varargin{2};
    V3 = varargin{3};
    A1 = doublearea([ P;V2;V3],[1:n;n+[1:n;n+(1:n)]]');
    A2 = doublearea([V1; P;V3],[1:n;n+[1:n;n+(1:n)]]');
    A3 = doublearea([V1;V2; P],[1:n;n+[1:n;n+(1:n)]]');
    A  = doublearea([V1;V2;V3],[1:n;n+[1:n;n+(1:n)]]');
    if size(P,2)>2 && max(abs(sum([A1 A2 A3],2)-A))>1e-14
      warning('Possibly negative coordinates. Not supported in dim~=2');
    end
    B = bsxfun(@rdivide,[A1 A2 A3],A);
  otherwise
    error(sprintf('%d-simplices not supported',numel(varargin)));
  end
end
