function rescaleOFF2(original_filename, rescaled_filename)
  % Takes an OFF file and rescales it so all vertex positions are between 0
  % and 1 and writes it to a new file

 mesh = readOFF(original_filename); 
 
 % center
 B = sum(mesh.V)/size(mesh.V,1);
 mesh.V = mesh.V - repmat(B,size(mesh.V,1),1);
  
 % translate to origin
 mesh.V = mesh.V-min(min(mesh.V));
 % scale to fit in unit box
 mesh.V = mesh.V./max(max(mesh.V));

 writeOFF(rescaled_filename,mesh);
end
