function [V,T,F] = readMESH( filename )
  % readMESH reads an MESH file with vertex/face/tet information
  %
  % [V,T,F] = readMESH( filename )
  %
  % Input:
  %  filename  path to .mesh file
  % Outputs:
  %  V  #V by 3 list of vertices
  %  T  #T by 4 list of tet indices
  %  F  #F by 3 list of triangle indices
  %

  fp = fopen(filename,'r');

  % First line is mandatory header
  MESHheader = eat_comments(fp,'#');
  MESHheader = regexprep(MESHheader,'[\n\r]+','');
  if(strcmp(MESHheader,'MeshVersionFormatted 1')==0)
    error('First line should be "MeshVersionFormatted 1" not ("%s")...',MESHheader);
  end
  % force read line feed
  fscanf(fp,'\n');
  Dimension3 = eat_comments(fp,'#');
  Dimension3 = regexprep(Dimension3,'[\n\r]+','');
  % second line is mandatory Dimension 3
  if(strcmp(Dimension3,'Dimension 3')==0)
    % tetgen likes to put the 3 on the next line
    % try to append next word hoping its a 3
    % force read line feed
    Dimension3 = [Dimension3 ' ' eat_comments(fp,'#')];
    if(strcmp(Dimension3,'Dimension 3')==0)
      error('Second line should be "Dimension 3"...');
    end
  end
  Vertices = eat_comments(fp,'#');
  Vertices = regexprep(Vertices,'[\n\r]+','');
  % thrid line is mandatory Vertices
  if(strcmp(Vertices,'Vertices')==0)
    error('Third line should be "Vertices"...');
  end
  
  % read vertex count
  if ispc % is windows
      num_vertices = fscanf(fp,'%d\n\r',1);
  else
      num_vertices = fscanf(fp,'%d\n',1);
  end
  % read num_vertices many sets of vertex coordinates (x,y,z,ref)
  V = fscanf(fp,'%g',4*num_vertices);
  V = reshape(V,4,num_vertices)';
  V = V(:,1:3);

  Triangles = eat_comments(fp,'#');
  Triangles = regexprep(Triangles,'[\n\r]+','');
  if(strcmp(Triangles,'Edges')==1)
    % read triangle count
    if ispc % is windows
        num_edges = fscanf(fp,'%d\n\r',1);
    else
        num_edges = fscanf(fp,'%d\n',1);
    end
    % read num_triangles many sets of face indices (a,b,c,ref)
    E = fscanf(fp,'%d',3*num_edges);
    E = reshape(E,3,num_edges)';
    E = E(:,1:2);
    Triangles = eat_comments(fp,'#');
    Triangles = regexprep(Triangles,'[\n\r]+','');
  end


  % forth non numeric is mandatory Triangles
  if(strcmp(Triangles,'Triangles')==0)
    error('Fourth (non-number) line should be "Triangles"...');
  end
  % read triangle count
  if ispc % is windows
      num_triangles = fscanf(fp,'%d\n\r',1);
  else
      num_triangles = fscanf(fp,'%d\n',1);
  end
  % read num_triangles many sets of face indices (a,b,c,ref)
  F = fscanf(fp,'%d',4*num_triangles);
  F = reshape(F,4,num_triangles)';
  F = F(:,1:3);

  Tetrahedra = eat_comments(fp,'#');
  Tetrahedra = regexprep(Tetrahedra,'[\n\r]+','');
  % forth non numeric is mandatory Tetrahedra
  if(strcmp(Tetrahedra,'Tetrahedra')==0)
    error('Fifth (non-number) line should be "Tetrahedra"...');
  end
  % read tetrahedra count
  
  if ispc % is windows
      num_tetrahedra = fscanf(fp,'%d\n\r',1);
  else
      num_tetrahedra = fscanf(fp,'%d\n',1);
  end
  % read num_tetrahedra many sets of tet indices (a,b,c,d,ref)
  T = fscanf(fp,'%d',5*num_tetrahedra);
  T = reshape(T,5,num_tetrahedra)';
  T = T(:,1:4);
  
  fclose(fp);
end
