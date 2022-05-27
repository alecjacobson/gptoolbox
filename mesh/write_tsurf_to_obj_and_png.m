function write_tsurf_to_obj_and_png(fileprefix, tsh)

  % Check that plot is set up for per-vertex data
  assert(strcmp(tsh.FaceColor,'interp'));
  % Even though you may pass it  with 'CData' seems it gets moved to
  % 'FaceVertexCData'
  assert(size(tsh.FaceVertexCData,2) == 1);
  assert(size(tsh.FaceVertexCData,1) == size(tsh.Vertices,1));

  % get color map
  CM = tsh.Parent.Colormap;
  im = permute(flipud(CM),[1 3 2]);
  % make it a square
  im = repmat(im,[1 size(im,1), 1]);
  % use at least 1024
  im = imresize(im,ceil(1024/size(im,1)),'nearest');
  imwrite(im,[fileprefix '.png']);

  V = tsh.Vertices;
  F = tsh.Faces;
  CA = caxis(tsh.Parent);
  D = min( max( (tsh.FaceVertexCData-CA(1))/(CA(2)-CA(1)) , 0+1/size(im,1)),1-1/size(im,1));

  save([fileprefix '.mat'],'V','F','CA','D','CM');

  writeOBJ([fileprefix '.obj'],V,F,[D D],F);
  text = fileread([fileprefix '.obj']);
  fp = fopen([fileprefix '.obj'],'w');
  fprintf(fp,'mtllib %s.mtl\n',fileprefix);
  fprintf(fp,'usemtl Textured\n');
  fprintf(fp,'%s',text);
  fclose(fp);

  fp = fopen([fileprefix '.mtl'],'w');
  fprintf(fp,'newmtl Textured\n');
  fprintf(fp,'Ka 1.000 1.000 1.000\n');
  fprintf(fp,'Kd 1.000 1.000 1.000\n');
  fprintf(fp,'Ks 0.000 0.000 0.000\n');
  fprintf(fp,'d 1.0\n');
  fprintf(fp,'illum 2\n');
  fprintf(fp,'# the diffuse texture map\n');
  fprintf(fp,'map_Kd %s.png\n',fileprefix);
  fclose(fp);


end
