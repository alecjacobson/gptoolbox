function [V,T,F] = regular_tetrahedral_mesh(varargin)
  % REGULAR_TETRAHEDRAL_MESH Generates a regular tetrahedral mesh with
  % dimensions (nx,ny,nz)
  %
  % [V,T,F] = regular_tetrahedral_mesh(nx,ny,nz)
  % [V,T,F] = regular_tetrahedral_mesh(nx)
  % [V,T,F] = regular_tetrahedral_mesh([nx,ny,nz])
  %
  % Input:
  %   nx  number of points in x direction on grid
  %   ny  number of points in y direction on grid
  %   nz  number of points in z direction on grid
  % Output:
  %   V  list of vertex coordinates in 3D
  %   T  tetrahedra list of indices into V
  %   F  triangle list of face indices into V
  % %   N  sparse adjacency matrix, #vertices by #vertices
  %
  % Example
  %    [V,T,F] = regular_tetrahedral_mesh(3,3,3);
  %    tetramesh(T,V); %also try tetramesh(T,V,'FaceAlpha',0);
  %    trisurf(F,V(:,1),V(:,2),V(:,3));
  % 
  % See also delaunayn, tetramesh
  %
  if nargin==1
    if numel(varargin{1})==1
      nx = varargin{1};
      ny = nx;
      nz = nx;
    else
      nx = varargin{1}(1,1);
      ny = varargin{1}(1,2);
      nz = varargin{1}(1,3);
    end
  else
    nx = varargin{1};
    ny = varargin{2};
    nz = varargin{3};
  end

  % Create a grid
  [x,y,z] = meshgrid(linspace(0,1,nx),linspace(0,1,ny),linspace(0,1,nz));
  V = [x(:) y(:) z(:)];
  % meshgrid flips x and y ordering
  idx = reshape(1:prod([ny,nx,nz]),[ny,nx,nz]);
  v1 = idx(1:end-1,1:end-1,1:end-1);v1=v1(:);
  v2 = idx(1:end-1,2:end,1:end-1);v2=v2(:);
  v3 = idx(2:end,1:end-1,1:end-1);v3=v3(:);
  v4 = idx(2:end,2:end,1:end-1);v4=v4(:);
  v5 = idx(1:end-1,1:end-1,2:end);v5=v5(:);
  v6 = idx(1:end-1,2:end,2:end);v6=v6(:);
  v7 = idx(2:end,1:end-1,2:end);v7=v7(:);
  v8 = idx(2:end,2:end,2:end);v8=v8(:);
  T = [ ...
    v1  v3  v8  v7; ...
    v1  v8  v5  v7; ...
    v1  v3  v4  v8; ...
    v1  v4  v2  v8; ...
    v1  v6  v5  v8; ...
    v1  v2  v6  v8];

  % delaunayn IS BROKEN, REALLY SHOULD NOT USE THIS
  %if(nx==2 && ny==2 && nz==2)
  %  warning(['delaunayn vomits on 2x2x2...' ...
  %    'adding center point at [0.5,0.5,0.5].']);
  %  V = [V;0.5,0.5,0.5];
  %end
  %T = delaunayn(V);

  % Much faster (~17x) than calling boundary_faces
  if nargout>2
    I = bsxfun(@plus,1:ny,(0:nx*ny:(nx*ny*(nz-1)))');
    ISE = I(1:end-1,1:end-1);
    INE = I(1:end-1,2:end);
    ISW = I(2:end,1:end-1);
    INW = I(2:end,2:end);
    F1 = [ISE(:) INW(:) ISW(:);ISE(:) INE(:) INW(:)];
    F1 = [F1;fliplr(F1)+ny*(nx-1)];

    I = bsxfun(@plus,1:ny:(nx-1)*ny+1,(0:nx*ny:(nx*ny*(nz-1)))');
    ISE = I(1:end-1,1:end-1);
    INE = I(1:end-1,2:end);
    ISW = I(2:end,1:end-1);
    INW = I(2:end,2:end);
    F2 = [ISE(:) INW(:) ISW(:);ISE(:) INE(:) INW(:)];
    F2 = [F2;fliplr(F2)+(ny-1)];

    I = bsxfun(@plus,1:ny:(nx-1)*ny+1,(0:ny-1)');
    ISE = I(1:end-1,1:end-1);
    INE = I(1:end-1,2:end);
    ISW = I(2:end,1:end-1);
    INW = I(2:end,2:end);
    F3 = [ISE(:) INW(:) ISW(:);ISE(:) INE(:) INW(:)];
    F3 = [F3;fliplr(F3)+nx*ny*(nz-1)];

    F = [F1;fliplr(F2);F3];
    %BF = boundary_faces(T);
    %setdiff(sort(F,2),sort(BF,2),'rows')
    %setdiff(sort(BF,2),sort(F,2),'rows')
  end

  % determine neighbors using fact that vertices are in order
%   I = (1:size(V,1))';
%   J = [ ...
%     I+1 I-1 ...
%     I+nx*ny I-nx*ny ...
%     I+nx*ny+1 I-nx*ny-1 ...
%     I+ny I-ny ...
%     I+ny+1 I-ny-1 ...
%     I+ny+nx*ny I-ny-nx*ny ...
%     I+ny+1+nx*ny I-ny-1-nx*ny];
%   C = [ ...
%     mod(I,4) > 0  mod(I,4) ~= 1 ...
%     I<=(nx*ny*nz-nx*ny) I>(nx*ny) ...
%     (I<=(nx*ny*nz-nx*ny) & mod(I,4) > 0) (I>(nx*ny) & mod(I,4)~= 1) ...
%     (mod(I-1,nx*ny) < nx*ny-ny) (mod(I-1,nx*ny) >= ny) ...
%     (mod(I-1,nx*ny) < nx*ny-ny & mod(I,4) > 0) ...
%     (mod(I-1,nx*ny) >= ny & mod(I,4)~= 1) ...
%     (I<=(nx*ny*nz-nx*ny) & mod(I-1,nx*ny) < nx*ny-ny)  ...
%     (I>(nx*ny) & mod(I-1,nx*ny) >= ny) ...
%     (I<=(nx*ny*nz-nx*ny) & mod(I-1,nx*ny) < nx*ny-ny & mod(I,4) > 0) ...
%     (I>(nx*ny) & mod(I-1,nx*ny) >= ny & mod(I,4) ~= 1)];
%   I = repmat(I,1,14);
  
  %N = sparse([I(V) J(V)],[J(V) I(V)],[V(V) V(V)]*0.5,size(V,1),size(V,1));
  

  % replace 

  % This is a bad idea because then memory layout doesn't match regular grid
  % (so finite difference won't match)
  %% convenient to have these reordered so that vertices that occur only on the
  %% surface come first
  %[V,T,F] = faces_first(V,T,F);
  %%N = sparse([IM(I(C)) IM(J(C))],[IM(J(C)) IM(I(C))],[C(C) C(C)]*0.5,size(C,1),size(C,1));
  
end

