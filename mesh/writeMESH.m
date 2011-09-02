function writeMESH( filename, V,T,F)
  % WRITEMESH writes an MESH file with vertex/face/tet information
  %
  % writeMESH(filename,V,T,F)
  %
  %
  % Input:
  %  filename  path to .mesh file
  %  V  #V by 3 list of vertices
  %  T  #T by 4 list of tet indices
  %  F  #F by 3 list of triangle indices
  %
  % See also readMESH
  %
  disp(['writing: ',filename]);
  fp = fopen(filename,'w');
  % mandatory header info
  fprintf(fp,'MeshVersionFormatted 1\n');
  fprintf(fp,'Dimension 3\n');
  % vertices header
  fprintf(fp,'Vertices\n');
  % number of vertices
  fprintf(fp,'%d\n',size(V,1));
  % vertex positions
  fprintf(fp,'%0.15g %0.15g %0.15g 1\n',V');
  % triangles header
  fprintf(fp,'Triangles\n');
  % number of triangles
  fprintf(fp,'%d\n',size(F,1));
  % triangle indices
  fprintf(fp,'%d %d %d 1\n',F');
  % tetrahedra header
  fprintf(fp,'Tetrahedra\n');
  % number of tetrahedra 
  fprintf(fp,'%d\n',size(T,1));
  % tetrahedra indices
  fprintf(fp,'%d %d %d %d 1\n',T');
  % end
  fprintf(fp,'End');
  fclose(fp);
end

