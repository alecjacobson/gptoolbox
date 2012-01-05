function h = animated_tetramesh(T,X,Y,Z)
  % ANIMATED_TETRAMESH animated, interactive tetramesh plot
  %
  % h = animated_tetramesh(F,X,Y,Z)
  %
  % Inputs:
  %   F  #F by 4 list of tetrahedra indices
  %   V  #V by dim by #{frames | 1} list of vertex positions
  %    or
  %   X  #V by #{frames | 1} list of X positions
  %   Y  #V by #{frames | 1} list of Y positions
  %   Z  #V by #{frames | 1} list of Z positions
  % Outputs:
  %   h  handle to plot
  %
  % See also: tsurf, trisurf, tetramesh, animated_trisurf
  %

  % this is not working because tetramesh just displays #T surfaces
  % I guess I could do that: just build a "surface" mesh with faces for each
  % face in the tet mesh, then I could tweak Vertices data on the fly, in fact
  % then its just a wrapper for animated_trisurf
  % get all faces

  allF = [ ...
    T(:,1) T(:,2) T(:,3); ...
    T(:,1) T(:,3) T(:,4); ...
    T(:,1) T(:,4) T(:,2); ...
    T(:,2) T(:,4) T(:,3)];
  % sort rows so that faces are reorder in ascending order of indices
  sortedF = sort(allF,2);
  % determine uniqueness of faces
  [F,m,n] = unique(sortedF,'rows');

  if nargin == 2
    V = X;
    dim = size(V,2);
    X = squeeze(V(:,1,:));
    Y = squeeze(V(:,2,:));
    if dim == 3
      Z = squeeze(V(:,3,:));
    else
      Z = 0*X(:,1);
    end
  end

  h = animated_trisurf(F,X,Y,Z);

end

