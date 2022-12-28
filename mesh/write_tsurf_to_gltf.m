function write_tsurf_to_gltf(filename, tsh)
  drawnow;
  % Check that plot is set up for per-vertex data
  assert(strcmp(tsh.FaceColor,'interp'));
  % Even though you may pass it  with 'CData' seems it gets moved to
  % 'FaceVertexCData'
  assert(size(tsh.FaceVertexCData,2) == 1);
  assert(size(tsh.FaceVertexCData,1) == size(tsh.Vertices,1));

  % get color map
  CM = tsh.Parent.Colormap;
  CA = caxis(tsh.Parent);
  % pad colormap to power of 2
  pad = 2^ceil(log2(size(CM,1))) - size(CM,1);
  CA = [CA(1) CA(1) + (CA(2)-CA(1)) * (size(CM,1)+pad)/size(CM,1)];
  CM(end+(1:pad),:) = repmat(CM(end,:),pad,1);
  im = permute(flipud(CM),[1 3 2]);

  D = (tsh.FaceVertexCData-CA(1))/(CA(2)-CA(1));
  R = axisangle2matrix([1 0 0],pi/2);
  VN = -tsh.VertexNormals;
  if ~isempty(VN);
   VN = VN*R;
  end

  writeGLTF_and_validate( ...
    filename, ...
    tsh.Vertices * R , ...
    tsh.Faces, ...
    'Normals',VN, ...
    'MagFilter','NEAREST', ...
    'MinFilter','NEAREST', ...
    'TextureCoordinates',[D D], ... 
    'TextureImage',im);
end

