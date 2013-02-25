function [V,F,UV,TF,N,NF] = readOBJ( filename )
  % READOBJ reads an OBJ file with vertex/face information
  %
  % [V,F,UV,TF,N,NF] = readOBJ( filename )
  %
  % Input:
  %  filename  path to .obj file
  % Outputs:
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %  UV  #V by 2 list of texture coordinates
  %  TF  #F by 3 list of triangle texture coordinates
  %  N  #V by 3 list of normals
  %  NF  #F by 3 list of triangle corner normal indices into N
  %
  %
  % WARNING: This is at least 40 times slower than readOFF but probably much much
  % slower... Because it's probably quadratic
  %
  % See also: load_mesh, readOBJfast, readOFF




V = [];
UV = [];
F = [];
TF = [];
N = [];
NF = [];
triangulated = false;
fp = fopen( filename, 'r' );
type = fscanf( fp, '%s', 1 );
count = 0;
while strcmp( type, '' ) == 0
    line = fgets(fp);
    if strcmp( type, 'v' ) == 1
        v = sscanf( line, '%lf %lf %lf' );
        V = [V; v'];
    elseif strcmp( type, 'vt')
        v = sscanf( line, '%f %f %f' );
        UV = [UV; v'];
    elseif strcmp( type, 'vn')
        n = sscanf( line, '%f %f %f' );
        N = [N; n'];
    elseif strcmp( type, 'f' ) == 1
        [t, count] = sscanf(line, '%d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d');
        if (count>2)
            tf = t(2:3:end);
            nf = t(3:3:end);
            t = t(1:3:end);
        else
            [t, count] = sscanf(line, '%d/%d %d/%d %d/%d %d/%d %d/%d');
            if (count>2)
                tf = t(2:2:end);
                t = t(1:2:end);
                nf = -ones(numel(t),1);
            else
              [t, count] = sscanf(line, '%d//%d %d//%d %d//%d %d//%d %d//%d');
              if (count>2)
                  nf = t(2:2:end);
                  t = t(1:2:end);
                  tf = -ones(numel(t),1);
              else
                  [t, count] = sscanf( line, '%d %d %d %d %d %d %d %d %d %d %d\n' );
                  tf = -ones(numel(t),1);
                  nf = -ones(numel(t),1);
              end
            end
        end
        assert(numel(t) == numel(tf));
        assert(numel(t) == numel(nf));
        if numel(t) > 3
          if ~triangulated
            warning('Trivially triangulating high degree facets');
          end
          triangulated = true;
        end
        for j = 2:numel(t)-1
          F = [F; t([1 j j+1])'];
          TF = [TF; tf([1 j j+1])'];
          NF = [NF; nf([1 j j+1])'];
        end
    elseif strcmp( type, '#' ) == 1
        %fscanf( fp, '%s\n', 1 );
        % ignore line
    end
    type = fscanf( fp, '%s', 1 );
end
fclose( fp );

%try
%    F = cell2mat(F);
%end

%% transform into array if all faces have the same number of vertices

if (size(UV,1)>0) UV = UV; end

end
