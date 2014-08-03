function writePLY(filename,V,F,mode);
  % WRITEPLY wrapper for write_ply
  %
  % Inputs:
  %   filename  path to output .ply file
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh triangle indices
  %   mode  followed by
  %      {'ascii'}                ASCII text data
  %      'binary_little_endian'   binary data, little endian
  %      'binary_big_endian'      binary data, big endian
  %
  %
  write_ply(V,F,filename,mode);
end
