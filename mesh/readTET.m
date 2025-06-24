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
  if numel(TETheader)>=3 && TETheader(end-2:end) == 'TET'
    d = fscanf(fp, '%d', 2);
    nV = d(1);
    nT = d(2);
    V = fscanf(fp,'%g %g %g',[3 nV])';
    T = fscanf(fp,'%d %d %d %d',[4 nT])'+1;
  else
    warning('no tet header, trying again without reading header') 
    % marco attene seems to use a different header in CDT
    fclose(fp);
    fp = fopen( filename, 'r' );
    % get line
    line = fgetl(fp);
    % read first number in line
    nV = sscanf(line, '%d', 1);
    % get line
    line = fgetl(fp);
    % read first number in line
    nT = sscanf(line, '%d', 1);
    V = fscanf(fp,'%g %g %g',[3 nV])';
    T = fscanf(fp,'%d %d %d %d %d',[5 nT])'+1;
    T = T(:,2:5);
  end


  fclose(fp);
  
  %disp('  - done.');
end

