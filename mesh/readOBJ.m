function [V,F,UV] = readOBJ( filename )
  % READOBJ reads an OBJ file with vertex/face information
  %
  % [V,F,UV] = readOBJ( filename )
  %
  % Input:
  %  filename  path to .obj file
  % Outputs:
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %  UV  #V by 2 list of texture coordinates
  %
  %
  % WARNING: This is at least 40 times slower than readOFF but probably much much
  % slower...
  %
  % See also: load_mesh, readOBJfast, readOFF




V = [];
UV = [];
F = [];
fp = fopen( filename, 'r' );
type = fscanf( fp, '%s', 1 );
while strcmp( type, '' ) == 0
    if strcmp( type, 'v' ) == 1
        v = fscanf( fp, '%g %g %g\n' );
        V = [V; v'];
    elseif strcmp( type, 'vt')
        v = fscanf( fp, '%g %g %g\n' );
        UV = [UV; v'];
    elseif strcmp( type, 'f' ) == 1
        line = fgets(fp);
        [t, count] = sscanf(line, '%d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d');

        if (count>2)
            t = t(1:3:end);
        else
            [t, count] = sscanf(line, '%d/%d %d/%d %d/%d %d/%d %d/%d');
            if (count>2)
                t = t(1:2:end);
            else
                [t, count] = sscanf( line, '%d %d %d %d %d %d %d %d %d %d %d\n' );
            end
        end
        F = [F; t'];
    elseif strcmp( type, '#' ) == 1
        fscanf( fp, '%s\n', 1 );
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
