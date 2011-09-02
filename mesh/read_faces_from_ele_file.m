function F = read_faces_from_ele_file(ele_file_name)
  % read face connectivity information from a standard .ele file (output of
  % executing 'triangle'
  %
  % Input:
  %   ele_file_name: name of .ele file
  %

  ele_file_handle = fopen(ele_file_name);
  header = fscanf(ele_file_handle,'%d %d %d',3);
  % read each face (prefixed by face index
  F = fscanf(ele_file_handle,'%d %d %d %d',[header(2)+1,header(1)])';
  % remove column of face indices
  F = F(:,2:end);
  fclose(ele_file_handle);
end
