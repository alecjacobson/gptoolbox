function [V,F,UV,C,N] = readOFF( filename )
  % READOFF reads an OFF file with vertex/face information
  %
  % [V,F,UV,C,N] = readOFF( filename )
  %
  % Input:
  %  filename  path to .obj file
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
  if (OFFheader(end-2:end) ~= 'OFF') warning('no OFF file!'); return; end
  OFFdim = 3;
  OFF_N = 0; OFF_C=0; OFF_ST=0;
  
  if find(OFFheader=='N') OFFdim = OFFdim+3; OFF_N=1; end
  if find(OFFheader=='C') OFFdim = OFFdim+3; OFF_C=1; end
  if find(OFFheader=='S') OFFdim = OFFdim+2; OFF_ST=1; end
  
  d = fscanf( fp, '%d', 3);
  nV = d(1); nF = d(2); nE = d(3);
  
  disp(sprintf('  - Reading %d vertices', nV));
  
  switch OFFdim
      case  3; OFFV = textscan( fp, '%f %f %f', nV);
      case  5; OFFV = textscan( fp, '%f %f %f %f %f', nV);
      case  6; OFFV = textscan( fp, '%f %f %f %f %f %f', nV);
      case  7; OFFV = textscan( fp, '%f %f %f %f %f %f %f', nV);
      case  8; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f', nV);
      case  9; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f %f', nV);
      case 10; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f %f %f', nV);
      case 11; OFFV = textscan( fp, '%f %f %f %f %f %f %f %f %f %f %f', nV);
      otherwise; error('Unsupported number of vertex entries');
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
    disp(sprintf('  - Reading %d faces', nF));
    F = cell2mat( textscan( fp, '%d %d %d %d %d %d', nF ) );
    F =  double(F(:, 2:(F(1,1)+1) ) + 1 );
  else
    F = [];
  end
  
  disp('  - done.');
end
