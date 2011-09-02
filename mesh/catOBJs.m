function catOBJs(input_files, output_file)
  [V,F] = readOBJ(input_files{1});
  for i = 2:size(input_files,1)
    [Vi,Fi] = readOBJ(input_files{i});
    [V, F] = cat_meshes(V,F,Vi,Fi);
  end
  writeOBJ(output_file,V,F);
end
