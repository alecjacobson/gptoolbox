function writeOBJ(filename, V,F,UV,N)
  % WRITEOBJ writes an OBJ file with vertex/face information
  %
  % writeOBJ(filename,V,F,UV,N)
  %
  % Input:
  %  filename  path to .obj file
  % Outputs:
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %  UV  #V by 2 list of texture coordinates
  %  N  #V by 3 list of normals
  %

disp(['writing: ',filename]);
f = fopen( filename, 'w' );

for k=1:size(V,1)
    fprintf( f, 'v %f %f %f\n', V(k,1), V(k,2), V(k,3) );
end

hasN =  exist('N','var') && ~isempty(N);
hasUV = exist('UV','var') && ~isempty(UV);

if hasUV
    for k=1:size(UV,1)
        fprintf( f, 'vt %f %f\n', UV(k,1), UV(k,2) );
    end
end

if hasN
    for k=1:size(V,1)
        fprintf( f, 'vn %f %f %f\n', N(k,1), N(k,2), N(k,3) );
    end
end
    
    
if ( (~hasN) && (~hasUV) )
    for k=1:size(F,1)
        fprintf( f, 'f %d %d %d\n', ...
            F(k,1), F(k,2), F(k,3));
    end
elseif ( (~hasN) && (hasUV) )
    for k=1:size(F,1)
        fprintf( f, 'f %d/%d %d/%d %d/%d\n', ...
            F(k,1), F(k,1), ...
            F(k,2), F(k,2), ...
            F(k,3), F(k,3) );
    end
elseif ( (hasN) && (~hasUV) )
    for k=1:size(F,1)
        fprintf( f, 'f %d//%d %d//%d %d//%d\n', ...
            F(k,1), F(k,1), ...
            F(k,2), F(k,2), ...
            F(k,3), F(k,3) );
    end
elseif ( (hasN) && (hasUV) )
    for k=1:size(F,1)
        fprintf( f, 'f %d/%d/%d %d/%d/%d %d/%d/%d\n', ...
            F(k,1), F(k,1), F(k,1),...
            F(k,2), F(k,2), F(k,2),...
            F(k,3), F(k,3), F(k,3) );
    end
end


fclose(f);
