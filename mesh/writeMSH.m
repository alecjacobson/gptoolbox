function writeMSH(filename,V,T,F)
  % WRITEMSH  write a tet mesh to gmsh's .msh file format
  %
  % writeMSH(filename,V,T,F)
  %
  % Inputs:
  %   V  #V by dim=3 list of vertex positions
  %   T  #T by 4 list of tet indices
  %   F  #F by 3 list of triangle indices

  %disp(['writing: ',filename]);
  fp = fopen(filename,'w');
  % Write header
  fprintf(fp,'$MeshFormat\n');
  % http://www.manpagez.com/info/gmsh/gmsh-2.2.2/gmsh_63.php
  fprintf(fp,'%d %d %d\n',2,0,8);
  fprintf(fp,'$EndMeshFormat\n');
  % write nodes
  fprintf(fp,'$Nodes\n');
  fprintf(fp,'%d\n',size(V,1));
  fprintf(fp,'%d %0.17g %0.17g %0.17g\n',[1:size(V,1);V']);
  fprintf(fp,'$EndNodes\n');
  % write elements
  fprintf(fp,'$Elements\n');
  fprintf(fp,'%d\n',size(T,1)+size(F,1));
  % write tets
  if ~isempty(T)
    fprintf(fp,'%d 4 2 0 1 %d %d %d %d\n',[1:size(T,1);T']);
  end
  if ~isempty(F)
    fprintf(fp,'%d 2 2 0 1 %d %d %d\n',[size(T,1)+(1:size(F,1));F']);
  end
  fprintf(fp,'$EndElements\n');
  fclose(fp);
end
