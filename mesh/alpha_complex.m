function [AF,DF] = alpha_complex(varargin)
  % ALPHA_COMPLEX Compute the alpha complex for a given set of points and a
  % given alpha. That is, elements of the delaunay mesh for these point whose
  % circumsphere radius is at most 1/alpha
  %
  % [AF,DF] = alpha_complex(V,alpha)
  % [AF,DF] = alpha_complex(V,alpha,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by dim list of point positions
  %   alpha  aka 1/max radius
  %   Optional:
  %     'OnlyShrinkBoundary' followed by true or false, whether to only remove
  %       elements on the boundary (see Holes), default is false
  %     'Holes' Specifies internal holes, only realy useful if
  %       ShrinkBoundaryOnly is true. Followed by a #H by dim list of hole
  %       positions, should be inside convex hull of V
  % Outputs:
  %   AF  #AF by dim+1 list of element indices in alpha complex
  %   DF  #DF by dim+1 list of element indices in delaunay complex
  %

  V = varargin{1};
  alpha = varargin{2};
  
  % defaults
  only_shrink_boundary = false;
  H = [];

  ii = 3;
  while ii <= nargin
  switch varargin{ii}
  case 'OnlyShrinkBoundary'
    if (ii+1)<=nargin && ~ischar(varargin{ii})
      ii = ii + 1;
      only_shrink_boundary = varargin{ii};
    else
      only_shrink_boundary = true;
    end
  case 'Holes'
    assert((ii+1)<=nargin);
    ii = ii + 1;
    H = varargin{ii};
  otherwise
    error(['Unsupported option: ' varargin{ii}]);
  end
  ii = ii + 1;
  end


  if size(V,2) == 3 && all(V(:,3) == 0)
    warning('Ignoring all zero Z-coordinate, treating as 2D');
    V = V(:,1:2);
  end

  % compute delaunay mesh
  dim = size(V,2);
  switch dim
  case 2
    DF = delaunay(V(:,1),V(:,2));
    % triangle side lengths
    l = [ ...
      sqrt(sum((V(DF(:,2),:)-V(DF(:,3),:)).^2,2)) ...
      sqrt(sum((V(DF(:,3),:)-V(DF(:,1),:)).^2,2)) ...
      sqrt(sum((V(DF(:,1),:)-V(DF(:,2),:)).^2,2)) ...
      ];
    % compute element circumradii
    % http://en.wikipedia.org/wiki/Circumscribed_circle#Circumscribed_circles_of_triangles
    dblA = doublearea_intrinsic(l);
    r = prod(l,2)./dblA;
  case 3
    DF = delaunay(V(:,1),V(:,2),V(:,3));
    error
  otherwise
    DF = delaunayn(V);
    error
  end

  % initialize with delaunay mesh
  AF = DF;
  % remove any elements containing hole positions
  if ~isempty(H)
    I = in_face(V,DF,H);
    AF = AF(~any(I,1),:);
    r = r(~any(I,1),:);
  end

  % get while loop running
  old_AF_size = size(AF,1)+1;
  while(old_AF_size ~= size(AF,1))
    old_AF_size = size(AF,1);
    % remove elements that are too big 
    keep = r<(1/alpha);
    if only_shrink_boundary
      % ... and on boundary
      B = on_boundary(AF);
      tsurf(AF(B,:),V);
      % get boundary faces
      keep = keep | ~B;
    end
    AF = AF(keep,:);
    r = r(keep);
  end
end
