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
  %
  % See also: load_mesh, readOBJfast, readOFF



numv = 0;
numf = 0;
numuv = 0;
numtf = 0;
numn = 0;
numnf = 0;

% Amortized array allocation
V = zeros(10000,3);
F = zeros(10000,3);
UV = zeros(10000,3);
TF = zeros(10000,3);
N = zeros(10000,3);
NF = zeros(10000,3);

triangulated = false;
fp = fopen( filename, 'r' );
type = fscanf( fp, '%s', 1 );
count = 0;
while strcmp( type, '' ) == 0
    line = fgets(fp);
    if strcmp( type, 'v' ) == 1
        v = sscanf( line, '%lf %lf %lf' );
        numv = numv+1;
        if(numv>size(V,1))
          V = cat(1,V,zeros(10000,3));
        end
        V(numv,:) = [v'];
    elseif strcmp( type, 'vt')
        v = sscanf( line, '%f %f %f' );
        numuv = numuv+1;
        if size(UV,2)>2 && length(v) == 2
            UV = UV(:,1:2);
        end
        if(numuv>size(UV,1))
          UV = cat(1,UV,zeros(10000,length(v)));
        end
        UV(numuv,:) = [v'];
    elseif strcmp( type, 'vn')
        n = sscanf( line, '%f %f %f' );
        numn = numn+1;
        if(numn>size(N,1))
          N = cat(1,N,zeros(10000,3));
        end
        N(numn,:) = [n'];
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
          numf = numf+1;
          if(numf>size(F,1))
            F = cat(1,F,zeros(10000,3));
          end
          F(numf,:) = [t([1 j j+1])'];
          numtf = numtf+1;
          if(numtf>size(TF,1))
            TF = cat(1,TF,zeros(10000,3));
          end
          TF(numtf,:) = [tf([1 j j+1])'];
          numnf = numnf+1;
          if(numnf>size(NF,1))
            NF = cat(1,NF,zeros(10000,3));
          end
          NF(numnf,:) = [nf([1 j j+1])'];
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
V = V(1:numv,:);
F = F(1:numf,:);
UV = UV(1:numuv,:);
TF = TF(1:numtf,:);
N = N(1:numn,:);
NF = NF(1:numnf,:);

%% transform into array if all faces have the same number of vertices

if (size(UV,1)>0) UV = UV; end

end
