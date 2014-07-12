function bundle  =  readBDL(filename)
% READBDL Read mesh(es) from a .bdl bundle file
%
% bundle  =  readBDL(filename)
% 
% Input:
%   filename  path to .bdl file
% Outputs:
%   bundle  struct containing mesh data

  function mesh = readBDLMeshChunk(fp)
  % read a mesh chunk from a file stream 
  %
  % Usage:
  %    mesh = readMeshChunk( fp)
  %
  % mesh fields: 
  % V  (required)  vertices
  % F              faces
  % col            colors
  % baseColor      a single color for the whole mesh
  % UV             uv coords
          
      magicID   = fread(fp, 1, 'uint32');
          chunkSize = fread(fp, 1,  'int64');  
          nV        = fread(fp, 1, 'uint32'); 
          nF        = fread(fp, 1, 'uint32');
          nVperFace = fread(fp, 1, 'uint32'); 
          nCperPos  = fread(fp, 1, 'uint32'); 
          nCperUV   = fread(fp, 1, 'uint32'); 
          nCperCol  = fread(fp, 1, 'uint32'); 
          nCperNor  = fread(fp, 1, 'uint32'); 
          strip     = fread(fp, 1, 'uchar');
          mesh.baseColor = fread(fp, 4, 'float32');
      
          mesh.V   = fread(fp, [nCperPos,nV],'float32');
      mesh.V = mesh.V';
          mesh.UV  = fread(fp, [nCperUV,nV], 'float32');
      mesh.UV = mesh.V';
          mesh.col = fread(fp, [nCperCol,nV],'float32');
      mesh.col = mesh.col';
          if nCperNor ~= 0 && nCperNor ~= 3 
            error('nCperNor has to be 0 or 3')
          end
          if nCperNor == 3
            mesh.nor = fread(fp, [nCperNor,nV],'float32');
        mesh.nor = mesh.nor';
          end
          mesh.F = fread(fp,  [nVperFace,nF],'uint32');      
          mesh.F = mesh.F+1; 
  end


    
  fp = fopen(filename,'r');
  if fp == -1
     error('could not open file')
  end

  chunk = 1;  
  while ~feof(fp)
     magicID = fread(fp,1,'uint32');
     if feof(fp)
        break
     end
     bck =  -sizeof('uint32');
     fseek(fp, bck, 'cof');  
     if magicID == 1
    bundle.meshes{chunk} = readBDLMeshChunk(fp);
    chunk = chunk+1;
     else
        chunkSize = fread(fp, 1, 'uint64');
        fseek(fp, chunkSize, 'cof');
     end
  end
