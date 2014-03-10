function catOBJs(input_files, output_file)
  % Concatenate meshes obj files into one mesh
  %
  % catOBJs(input_files, output_file) 
  %
  % Inputs:
  %   input_files  cell of paths to input files
  % Output:
  %   output_file  path to output OBJ
  %
  [V,F] = load_mesh(input_files{1});
  for i = 2:numel(input_files)
    [Vi,Fi] = load_mesh(input_files{i});
    [V, F] = cat_meshes(V,F,Vi,Fi);
  end
  writeOBJ(output_file,V,F);
end
