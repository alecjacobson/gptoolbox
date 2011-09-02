function V = read_vertices_from_node_file(node_file_name)
  % read vertex position information from a standard .node file (output of
  % executing 'triangle'
  %
  % Input:
  %   node_file_name: name of .node file
  %

  node_file_handle = fopen(node_file_name);
  header = fscanf(node_file_handle,'%d %d %d %d',4);
  % read in each vertex position (prefixed by index and postfixed by border
  % marker)
  V = fscanf(node_file_handle,'%d %f %f %d',[4,header(1)])';
  % remove vertex index column and boundary marker column
  V = V(:,2:3);
  fclose(node_file_handle);
end
