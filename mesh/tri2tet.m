function tri2tet(input_name,output_name,skeleton_name,samples_per_edge)
  % TRI2TET Script to make 3d tet mesh from 2d triangle mesh with optional
  % constraint to place tet vertices at skeleton vertices and samples along
  % skeleton bones.
  %
  % tri2tet(input_name,output_name,skeleton_name,samples_per_edge)
  %
  % Input:
  %  input_name  name of file contaning input vertices and faces of surface,
  %    should end in .off or .obj
  %  output_name  name of output file to be written with 3D tet mesh, should
  %    end in .mesh
  %  skeleton_name  optional name of .tgf containing skeleton vertices and
  %    edges
  %  samples_per_edge  optional sampling of "bones", edges of skeleton stored
  %    in skeleton_name
  %

  % read input mesh
  [V,F] = load_mesh(input_name);

  IV = [];
  % read skeleton vertices if file given
  if(exist('skeleton_name','var'))
    [IV,E] = readTGF(skeleton_name);
    if(exist('samples_per_edge','var'))
      IV = sample_edges(IV,E,samples_per_edge);
    end
  end

  % use tetgen to mesh interior with tetrahedra
  %[V,T,F] = tetgen(V,F,IV);
  [V,T,F] = tetgen([V;IV],F);
  % write tet mesh to file
  writeMESH(output_name,V,T,F);

end
