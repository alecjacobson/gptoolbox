function writeMESH( filename, V,T,F)
  % WRITEMESH writes an MESH file with vertex/face/tet information
  %
  % writeMESH(filename,V,T,F)
  %
  %
  % Input:
  %  filename  path to .mesh file
  %  V  #V by 3 list of vertices
  %  T  #T by 4|5 list of tet indices (additional column is color index)
  %  F  #F by 3|4 list of triangle indices (additional column is color index) 
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
  if size(F,2) == 4
    % triangle indices + color reference
    fprintf(fp,'%d %d %d %d\n',F');
  else
    % triangle indices
    fprintf(fp,'%d %d %d 1\n',F');
  end
  % tetrahedra header
  fprintf(fp,'Tetrahedra\n');
  % number of tetrahedra 
  fprintf(fp,'%d\n',size(T,1));
  if size(T,2) == 5 
    % tetrahedra indices
    fprintf(fp,'%d %d %d %d %d\n',T');
  else
    % tetrahedra indices
    fprintf(fp,'%d %d %d %d 1\n',T');
  end
  % end
  fprintf(fp,'End');
  fclose(fp);
end

