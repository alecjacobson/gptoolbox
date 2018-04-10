function [V,T] = readTET( filename )
  % READTET reads an TET file with vertex/tet information
  %
  % [V,T] = readTET( filename )
  %
  % Input:
  %  filename  path to .tet file
  % Outputs:
  %  V  #V by 3 list of vertices
  %  T  #T by 3 list of triangle indices
  %
  % See also: load_mesh, readOBJfast, readOBJ

  V = [];
  T = [];
  
  fp = fopen( filename, 'r' );
  TETheader = upper(fscanf( fp, '%s\n', 1 ));
  if TETheader(end-2:end) ~= 'TET'
    warning('no tet file!') 
    fclose(fp);
    return;
  end
  d = fscanf(fp, '%d', 2);
  nV = d(1);
  nT = d(2);

  V = fscanf(fp,'%g %g %g',[3 nV])';
  T = fscanf(fp,'%d %d %d %d',[4 nT])'+1;


  fclose(fp);
  
  %disp('  - done.');
end

