function writeMESH( filename, V,T,F,E)
  % WRITEMESH writes an MESH file with vertex/face/tet information
  %
  % writeMESH(filename,V,T,F)
  % writeMESH(filename,V,T,F,E)
  %
  % Input:
  %  filename  path to .mesh file
  %  V  #V by 3 list of vertices
  %  T  #T by 4|5 list of tet indices (additional column is color index)
  %  F  #F by 3|4 list of triangle indices (additional column is color index) 
  %  Optional:
  %    E  #E by 2|3 list of edge indices (additional column is color index) 
  %
  % See also readMESH
  %
  if ~exist('E','var')
    E = [];
  end
      
  %disp(['writing: ',filename]);
  fp = fopen(filename,'w');
  % mandatory header info
  fprintf(fp,'MeshVersionFormatted 1\n');
  fprintf(fp,'Dimension 3\n');
  % vertices header
  fprintf(fp,'Vertices\n');
  % number of vertices
  fprintf(fp,'%d\n',size(V,1));
  % vertex positions
  fprintf(fp,'%0.15g %0.15g %0.15g 0\n',V');
  % triangles header
  fprintf(fp,'Triangles\n');
  % number of triangles
  fprintf(fp,'%d\n',size(F,1));
  if size(F,2) == 4
    % triangle indices + color reference
    fprintf(fp,'%d %d %d %d\n',F');
  else
    % triangle indices
    fprintf(fp,'%d %d %d 0\n',F');
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
    fprintf(fp,'%d %d %d %d 0\n',T');
  end
  % tetrahedra header
  fprintf(fp,'Edges\n');
  % number of tetrahedra 
  fprintf(fp,'%d\n',size(E,1));
  if size(E,2) == 3 
    % edge indices
    fprintf(fp,'%d %d %d\n',E');
  else
    % edge indices
    fprintf(fp,'%d %d 0\n',E');
  end
  % end
  fprintf(fp,'End');
  fclose(fp);
end

