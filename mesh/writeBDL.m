function writeBDL(filename,V,F)
  % WRITEBDL Write a mesh to a .bdl file
  % 
  % writeBDL(filename,V,F)
  %
  % Inputs:
  %   filename  path to .bdl file
  %   V  #V by dim list of mesh vertices
  %   F  #F by 3 list of triangle indices

  function headerStart = writeBDLChunkHeader(fp, chunktype, nBytes)
  % writes chunk header into file stream
  %
  % Usage:
  %   headerStart = writeChunkHeader( fp, chunktype, nBytes )
    switch chunktype
        case 'MESH_CHUNK'
            magicID = 1;
        case 'MESH_GROUP'
            magicID = 2;
        otherwise
            magicID = 0;
    end
    
    headerStart = ftell(fp);
    fwrite(fp, magicID, 'uint32');
    fwrite(fp, nBytes, 'int64');
  end

  function bundle = writeBDLMeshChunk(fp, mesh)
    % write a mesh chunk to a file stream that can be processed by the external
    % viewer.
    %
    % Usage:
    %    writeMeshChunk( fp, mesh )
    %
    % mesh fields: 
    % V  (required)  vertices
    % F              faces
    % col            colors
    % baseColor      a single color for the whole mesh
    % UV             uv coords


    [nV, nCperPos] = size(mesh.V);

    hasFaces = isfield(mesh,'F');
    if ( hasFaces )
        [nF, nVperFace] = size(mesh.F);
    else
        nF = 0; nVperFace = mesh.nVperFace;
    end

    nCperCol = 0;
    if isfield(mesh,'col') nCperCol = size(mesh.col,2); end

    baseColor = [1,1,1,1];
    if isfield(mesh,'baseColor') baseColor = mesh.baseColor; end

    nCperUV = 0;
    if isfield(mesh,'UV') nCperUV = size(mesh.UV,2); end

    nCperNormal = 0;
    strip = 0;

    headerStart = writeBDLChunkHeader(fp, 'MESH_CHUNK', 0); % chunk size filled in later
    chunkStart = ftell(fp);

    fwrite(fp, nV, 'uint32');
    fwrite(fp, nF, 'uint32');
    fwrite(fp, nVperFace, 'uint32');
    fwrite(fp, nCperPos, 'uint32');
    fwrite(fp, nCperUV, 'uint32');
    fwrite(fp, nCperCol, 'uint32');
    fwrite(fp, nCperNormal, 'uint32');

    fwrite(fp, strip, 'uchar');
    fwrite(fp, baseColor, 'float32');

    fwrite(fp, mesh.V', 'float32');
    if isfield(mesh,'UV') fwrite(fp, mesh.UV', 'float32'); end
    if isfield(mesh,'col') fwrite(fp, mesh.col', 'float32'); end
    if isfield(mesh,'F') fwrite(fp, mesh.F'-1, 'uint32');  % faces are 1-based in matlab

    chunkSize = ftell(fp) - chunkStart;

    fseek(fp, headerStart, 'bof');
    writeBDLChunkHeader(fp, 'MESH_CHUNK', chunkSize); 
    fseek(fp,0,'eof');

  end

  mesh.V = V;
  mesh.F = F;
  fp = fopen(filename,'w');
  writeBDLMeshChunk(fp,mesh);
  fclose(fp);
end
