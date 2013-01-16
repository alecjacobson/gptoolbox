function writeOBJ(filename, V,F,UV,TF,N,NF)
  % WRITEOBJ writes an OBJ file with vertex/face information
  %
  % writeOBJ(filename,V,F,UV,N)
  %
  % Input:
  %  filename  path to .obj file
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %  UV  #UV by 2 list of texture coordinates
  %  TF  #TF by 3 list of corner texture indices into UV
  %  N  #N by 3 list of normals
  %  NF  #NF by 3 list of corner normal indices into N
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

if hasUV && (~exist('TF') || isempty(TF))
    TF = F;
end
if hasN && (~exist('NF') || isempty(NF))
    NF = F;
end

for k=1:size(F,1)
    if ( (~hasN) && (~hasUV) ) || (any(TF(k,:)<=0,2) && any(NF(k,:)<=0,2))
        fprintf( f, 'f %d %d %d\n', ...
            F(k,1), F(k,2), F(k,3));
    elseif ( hasUV && (~hasN || any(NF(k,:)<=0,2)))
        fprintf( f, 'f %d/%d %d/%d %d/%d\n', ...
            F(k,1), TF(k,1), ...
            F(k,2), TF(k,2), ...
            F(k,3), TF(k,3) );
    elseif ( (hasN) && (~hasUV || any(TF(k,:)<=0,2)))
        fprintf( f, 'f %d//%d %d//%d %d//%d\n', ...
            F(k,1), NF(k,1), ...
            F(k,2), NF(k,2), ...
            F(k,3), NF(k,3) );
    elseif ( (hasN) && (hasUV) )
        assert(all(NF(k,:)>0));
        assert(all(TF(k,:)>0));
        fprintf( f, 'f %d/%d/%d %d/%d/%d %d/%d/%d\n', ...
            F(k,1), TF(k,1), NF(k,1),...
            F(k,2), TF(k,2), NF(k,2),...
            F(k,3), TF(k,3), NF(k,3) );
    end
end


fclose(f);
