function [V,F,VC,FC] = image_mesh(filename,m)
  % IMAGE_MESH Given a path to .png file, create a mesh within the α>½ region of
  % the image with triangles placed economically so that per-vertex colors will
  % look reasonable.
  % 
  % Input:
  %   filename  path to .png file
  %   m  number of output triangles
  % Outputs:
  %   V  #V by 2 list of 2D vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   VC  #V by 3 list of vertex colors
  %   FC  #F by 3 list of face colors
  %
  % See also: bwmesh
  %

  [im,M,A] = imread(filename);
  if ~isempty(M)
      im = ind2rgb(im,M);
  end
  if isempty(A)
    A = ones(size(im,1),size(im,2));
  end
  [V,F] = bwmesh(A,'SmoothingIters',2,'Tol',0.2);
  %[V,F] = bwmesh(A);
  [X,Y] = meshgrid((1:size(im,2))-0.5,(size(im,1):-1:1)-0.5);
  
  for c = 1:3 
    V(:,2+c) = interp2(X,Y,im2double(im(:,:,c)),V(:,1),V(:,2));
  end

  % weird trick to get better quality elements. Based on Delaunay lifting idea.
  %V(:,end+1) = (V(:,1).^2+V(:,2).^2)*0.001;

  [V,F] = decimate_libigl(V,F,m,'Method','qslim');

  VC = V(:,3:5);
  V = V(:,1:2);

  if nargout>3
    BC = barycenter(V,F);
    FC = zeros(size(F,1),3);
    for c = 1:3
      FC(:,c) = interp2(X,Y,im2double(im(:,:,c)),BC(:,1),BC(:,2));
    end
  end
end

