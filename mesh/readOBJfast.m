function [V,F] = readOBJfast(filename)
  % readOBJfast
  %
  % reads an OBJ file quickly, but OBJ file should be formated *simply*
  %
  % Input:
  %  filename  path to .obj file
  % Outputs:
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %
  % filename should be a text file containing first vertex positions where each
  % line is:
  % v x y z
  % then a list of faces where each line is:
  % f i j k
  % and exactly that. Any lines before vertices, or after vertices but before
  % faces are ignored. Comments throughout vetex and face lines will surely
  % break this.
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % For complete(r) support use readOBJ
  % 
  % See also readOBJ
  %
  fp = fopen(filename);
  % read in vertices
  V = [];
  while(true)
    V = fscanf(fp,' v %g %g %g ',inf);
    if(prod(size(V)) > 0)
      break;
    else
      line = fgets(fp);
      if(prod(size(line)) == 0)
        fclose(fp);
        error('Bad format... Try readOBJ...');
      end
    end
  end
  V = reshape(V,3,size(V,1)/3)';
  % read in faces
  F = [];
  while(true)
    F = fscanf(fp,' f %d %d %d ',inf);
    if(prod(size(F)) > 0)
      break;
    else
      line = fgets(fp);
      if(prod(size(line)) == 0)
        fclose(fp);
        error('Bad format... Try readOBJ...');
      end
    end
  end
  F = reshape(F,3,size(F,1)/3)';
  fclose(fp);
end
