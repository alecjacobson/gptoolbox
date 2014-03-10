function writeSMF(filename,V,F)
  % WRITESMF Write a mesh (V,F) to a .smf file.
  % 
  % writeSMF(filename,V,F)
  %
  % Input:
  %  filename  path to .obj file
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %
  % See also: writeOBJ.m
  %
  f = fopen( filename, 'w' );
  fprintf( f, 'v %0.17g %0.17g %0.17g\n', V');
  fprintf( f, 'f %d %d %d\n', F');
  fclose(f);
end
