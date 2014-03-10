function clean_tetgen_mesh(input_file,output_file)
  % CLEAN_TETGEN_MESH TetGen writes ALL triangles to the .mesh file, that is
  % all faces of the tetrahedra in the 3D mesh.
  % clean_tetgen_mesh(input_file,output_file) reads the vertices and tetrahedra
  % in input_file (ignoring triangles) and determines the surfaces triangles
  % then write the vertices, surfaces triangles and tetrahedra to output_file
  %
  % clean_tetgen_mesh(input_file,output_file)
  %
  % Input:
  %   input_file  path to .mesh file containing vertices and tetrahedra
  %   output_file  path to .mesh file to be written
  % 
  [V,T,F] = readMESH(input_file);
  F = boundary_faces(T);
  % reverse orientation of triangles, tetgen seems to have backwards
  % orientation
  F = [F(:,3) F(:,2) F(:,1)];
  % rearrange vertices so that vertices on surface come before internal
  % vertices
  [V,T,F] = faces_first(V,T,F);
  writeMESH(output_file,V,T,F);
end
