function BW = readVOX(filename)
  % READVOX  Read a 3D binary image (e.g. resulting from voxelization) from a
  % .vox file. Ignores the color palette at the end of the file.
  %
  % Inputs:
  %   filename  path to .vox file
  % Outputs:
  %   BW  ny by nx by nz  3D binary image
  %
  % Examples:
  %   writeOBJ('temp.obj',V,F);
  %   !./poly2vox temp.obj temp.vox -h -v4
  %   BW = readVOX('temp.vox');
  %   % poly2vox uses an odd method to find voxel size
  %   r = max(((max(V)-min(V))))/(max(size(BW))-0.1);
  %   [X,Y,Z] = meshgrid( ...
  %     (0.5:size(BW,2)-0.5)*r+min(V(:,1)), ...
  %     (0.5:size(BW,1)-0.5)*r+min(V(:,2)), ...
  %     max(V(:,3))-(size(BW,3)-0.5:-1:0.5)*r);
  %   BC = [X(:) Y(:) Z(:)];
  %   [iV,iQ] = voxel_surface(BW,'Centers',BC);
  %   [iV,IM] = remove_unreferenced(iV,iQ);
  %   iQ = IM(iQ);
  %   clf;
  %   hold on;
  %   trisurf(iQ,iV(:,1),iV(:,2),iV(:,3),'FaceAlpha',0.8,'EdgeAlpha',0.8,'FaceColor','r');
  %   tsurf(F,V,'FaceAlpha',0.8,'EdgeAlpha',0.8,'FaceColor','g');
  %   hold off;
  %   axis equal;
  %   iF = [iQ(:,[1 2 3]);iQ(:,[1 3 4])];
  % 

  f = fopen(filename, 'r');
  header = fread(f,3,'uint32=>double')';
  BW = fread(f,prod(header),'uchar');
  fclose(f);

  BW = reshape(BW,fliplr(header));
  BW = flipud(BW);
  BW = permute(BW,[1 3 2]);
  % Flip and convert to logical
  BW = BW==0;
end
