% writes chunk header into file stream
%
% Usage:
%   headerStart = writeChunkHeader( fp, chunktype, nBytes )

function headerStart = writeChunkHeader(fp, chunktype, nBytes)

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
