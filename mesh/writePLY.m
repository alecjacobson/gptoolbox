function writePLY(filename, V, F, mode, VColor)
  % WRITEPLY wrapper for write_ply
  %
  % Inputs:
  %   filename  path to output .ply file
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh triangle indices
  %   mode  followed by
  %      'ascii'                ASCII text data
  %      {'binary_little_endian'}   binary data, little endian
  %      'binary_big_endian'      binary data, big endian
  %   VColor #V array optional parameter to set the vertex color
  %
  
  if nargin < 4
    mode = 'binary_little_endian';
  end

  if nargin < 5
      VColor = [];
  end
  
  fid = fopen(filename,'wt');
  fprintf(fid,'ply\nformat %s 1.0\n',mode);
  fprintf(fid,'element vertex %d\n',size(V,1));
  fprintf(fid,'property double x\n');
  fprintf(fid,'property double y\n');
  fprintf(fid,'property double z\n');
  
  concatV = V;
  if ~isempty(VColor)
    fprintf(fid,'property uchar red\n');
    fprintf(fid,'property uchar green\n');
    fprintf(fid,'property uchar blue\n');
    concatV = [V,VColor];
  end
  fprintf(fid,'element face %d\n',size(F,1));
  fprintf(fid,'property list int int vertex_indices\n');
  fprintf(fid,'end_header\n');
  FF = [size(F,2)*ones(size(F,1),1) F-1];
  switch mode
  case 'ascii'
    % do nothing
    printString = '%0.17f %0.17f %0.17f\n';
    if ~isempty(VColor)
        printString = '%0.17f %0.17f %0.17f %d %d %d\n';
    end
    fprintf(fid,printString,concatV');
    format = [repmat('%d ',1,size(FF,2)) '\n'];
    fprintf(fid,format,FF');
  case {'binary_little_endian','binary_big_endian'}
    fclose(fid);
    switch mode
    case 'binary_little_endian'
      fid = fopen(filename,'a','ieee-le');
    case 'binary_big_endian'
      fid = fopen(filename,'a','ieee-be');
    end
    fwrite(fid,concatV','double');
    fwrite(fid,FF','int');
  otherwise
    error('Unsupported format: %s',mode);
  end
  fclose(fid);

  %% OLD 200x slower way (possible to produce slightly smaller files and control
  %% over precision)
  %write_ply(V,F,filename,mode);
end
