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

%disp(['writing: ',filename]);
f = fopen( filename, 'w' );

assert(size(V,2) == 3);
fprintf( f, 'v %0.17g %0.17g %0.17g\n', V');

hasN =  exist('N','var') && ~isempty(N);
hasUV = exist('UV','var') && ~isempty(UV);

if hasUV
    switch size(UV,2)
    case 2
      fprintf( f, 'vt %0.17g %0.17g\n', UV');
    case 3
      fprintf( f, 'vt %0.17g %0.17g %0.17g\n', UV');
    end
end

if hasN
    %for k=1:size(N,1)
    %    fprintf( f, 'vn %f %f %f\n', N(k,1), N(k,2), N(k,3) );
    %end
    fprintf( f, 'vn %0.17g %0.17g %0.17g\n', N');
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
