function [XYZV,I,BC] = readSVOX(filename)
  % [XYZV,I,BC] = readSVOX(filename)
  %
  % Inputs:
  %   filename  path  to .svox file
  % Outputs:
  %   XYZV  #XYZV by 4 long list of raw sparse data
  %   I  y by x by z binary voxel occpuancy
  %   BC  y by x by z grid of voxel barycenters
  %   
  %
  % [xsize][ysize][zsize]
  f = fopen(filename, 'r');
  header = fread(f,4,'uint32=>double')';
  magic_numbers = fread(f,5,'single=>double')';

  XYZV = fread(f,4*header(4),'int32=>double')';
  XYZV = reshape(XYZV+1,4,header(4))';
  sz = header([2 1 3]);
  I = zeros(sz);
  I(sub2ind(size(I),XYZV(:,2),XYZV(:,1),XYZV(:,3))) = XYZV(:,4);
  min_corner = [nan nan nan];
  h = nan;
  fclose(f);

  fx = magic_numbers(1);
  fy = magic_numbers(2);
  fz = magic_numbers(3);
  f = magic_numbers(4);
  fudge = magic_numbers(5);
  [X,Y,Z] = meshgrid( ...
     (0.5:size(I,2)-0.5)/f+fy,...
     (0.5:size(I,1)-0.5)/f+fx,...
     (0.5:size(I,3)-0.5)/f-fz);
  BC = [X(:) Y(:) Z(:)];

end
