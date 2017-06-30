function writeWIRE(filename, V,E)
  % WRITEWIRE writes an .wire file with vertex/edge information (see PyMesh
  % documentation).
  %
  % writeWIRE(filename,V,F,UV,N)
  %
  % Input:
  %  filename  path to .wire file
  %  V  #V by 3 list of vertices
  %  E  #E by 2 list of edge indices into E
  %

  f = fopen( filename, 'w' );
  if size(V,2) == 2
    warning('Appending 0s as z-coordinate');
    V(:,end+1:3) = 0;
  else
    assert(size(V,2) == 3);
  end

  fprintf( f, 'v %0.17g %0.17g %0.17g\n', V');
  assert(size(E,2) == 2,'Edge list should have two columns');
  fprintf( f, 'l %d %d\n', E');

  fclose(f);
end
