function bundle  =  loadBundle(filename)
% bundle  =  loadBundle(filename)
% read a bundle from a file 
		
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
