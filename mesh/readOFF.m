function [V,F,UV,C,N] = readOFF( filename )
  % READOFF reads an OFF file with vertex/face information
  %
  % [V,F,UV,C,N] = readOFF( filename )
  %
  % Input:
  %  filename  path to .off file
  % Outputs:
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %  UV  #V by 2 list of texture coordinates
  %  C  #V by 3 list of colors
  %  N  #V by 3 list of normals
  %
  % See also: load_mesh, readOBJfast, readOBJ

% (C) 2007 Denis Kovacs, NYU
%-------------------------------------------------------------------------

  V = [];
  F = [];
  UV = [];
  C = [];
  N = [];
  
  fp = fopen( filename, 'r' );
  OFFheader = upper(fscanf( fp, '%s\n', 1 ));
  if OFFheader(end-2:end) ~= 'OFF'
    warning('no OFF file!') 
    fclose(fp);
    return;
  end
  OFFdim = 3;
  OFF_N = 0; OFF_C=0; OFF_ST=0;
  
  if find(OFFheader=='N') OFFdim = OFFdim+3; OFF_N=1; end
  if find(OFFheader=='C') OFFdim = OFFdim+3; OFF_C=1; end
  if find(OFFheader=='S') OFFdim = OFFdim+2; OFF_ST=1; end

  % eat any comments before
  line = eat_comments(fp,'#');
  
  d = sscanf( line, '%d', 3);
  nV = d(1); nF = d(2); nE = d(3);
  
  %disp(sprintf('  - Reading %d vertices', nV));
  
  switch OFFdim
      case  3; OFFV = textscan( fp, '%f %f %f', nV);
      case  5; OFFV = textscan( fp, '%f %f %f %f %f', nV);
      case  6; OFFV = textscan( fp, '%f %f %f %f %f %f', nV);
      case  7; OFFV = textscan( fp, '%f %f %f %f %f %f %f', nV);
      case  8; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f', nV);
      case  9; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f %f', nV);
      case 10; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f %f %f', nV);
      case 11; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f %f %f %f', nV);
      otherwise; fclose(fp); error('Unsupported number of vertex entries');
  end
  
  try
     OFFV = cell2mat(OFFV); 
  end
  
  OFFdim = 1;
  V = OFFV(:,OFFdim:(OFFdim+2)); OFFdim = OFFdim + 3;
  if (OFF_N) N = OFFV(:,OFFdim:(OFFdim+2)); OFFdim = OFFdim + 3; end
  if (OFF_C) C = OFFV(:,OFFdim:(OFFdim+2)); OFFdim = OFFdim + 3; end
  if (OFF_ST) UV = OFFV(:,OFFdim:(OFFdim+1)); OFFdim = OFFdim + 2; end
  
  if (nF ~= 0)
    %disp(sprintf('  - Reading %d faces', nF));
    temp = textscan( fp, '%d %d %d %d %d %d %d %d %d %d %d', nF );
    sz = temp{1}(1);
    if all(sz == cell2mat(temp(1)))
      F = double (cell2mat( temp(2:sz+1 ))) +1;
    else
      warning('Trivially triangulating high degree facets');
      F = zeros(sum(temp{1}-2),3);
      fi = 1;
      for f = 1:size(temp{1},1)
        for j = 3:temp{1}(f);
          F(fi,:) = [temp{2}(f) temp{j}(f) temp{j+1}(f)]+1;
          fi = fi+1;
        end
      end
    end
  else
    F = [];
  end

  fclose(fp);
  
  %disp('  - done.');
end
