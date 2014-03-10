function writeOFF(filename, V,F,UV,C,N)
  % WRITEOFF writes an OFF file with vertex/face information
  %
  % writeOFF(filename,V,F)
  % writeOFF(filename,V,F,UV,C,N)
  %
  % Input:
  %  filename  path to .obj file
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %  UV  #V by 2 list of texture coordinates
  %  C  #V by 3 list of colors
  %  N  #V by 3 list of normals

%disp(['writing: ',filename]);

hasN =  exist('N','var') && ~isempty(N);
hasUV = exist('UV','var') && ~isempty(UV);
hasC = exist('C','var') && ~isempty(C);

OFFheader = 'OFF';

OFFV = V;

if hasN
    OFFheader = ['N',OFFheader];
    OFFV = [OFFV N];
end

if hasC
    OFFheader = ['C',OFFheader];
    OFFV = [OFFV C];
end

if hasUV
    OFFheader = ['ST',OFFheader];
    OFFV = [OFFV UV];
end 

f = fopen( filename, 'wt' );
fprintf(f, [OFFheader,'\n']);
fprintf(f, '%d %d 0\n', size(V,1), size(F, 1));

switch size(OFFV, 2)
    case  3; fprintf(f, '%0.17g %0.17g %0.17g\n', OFFV');
    case  5; fprintf(f, '%0.17g %0.17g %0.17g %0.17g %0.17g\n', OFFV');
    case  6; fprintf(f, '%0.17g %0.17g %0.17g %0.17g %0.17g %0.17g\n', OFFV');
    case  7; fprintf(f, '%0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g\n', OFFV');
    case  8; fprintf(f, '%0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g\n', OFFV');
    case  9; fprintf(f, '%0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g\n', OFFV');
    case 10; fprintf(f, '%0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g\n', OFFV');
    case 11; fprintf(f, '%0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g %0.17g\n', OFFV');
    otherwise; error('Unsupported number of vertex entries');
end

if (~isempty(F)) F; F = F - 1; end

switch size(F, 2)
    case 0;
    case 1; fprintf( f, '1 %d\n', F');
    case 2; fprintf( f, '2 %d %d\n', F');
    case 3; fprintf( f, '3 %d %d %d\n', F');
    case 4; fprintf( f, '4 %d %d %d %d\n', F');
    case 5; fprintf( f, '5 %d %d %d %d %d\n', F');
    otherwise; error('Unsupported number of vertex entries');
end    

fclose(f);
