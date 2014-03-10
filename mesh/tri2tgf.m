function tri2tgf(input_name,output_name)
  % TRI2TGF Script to convert a triangle mesh to a tgf graph (skeleton file)
  %
  % tri2tgf(input_name,output_name)
  %
  % Inputs:
  %  input_name  name of file contaning input vertices and faces of surface,
  %    should end in .off or .obj
  %  output_name  name of output file to be written with 3D tet mesh, should
  %    end in .mesh

  % read input mesh
  [V,F] = load_mesh(input_name);
  % Find all edges in mesh, note internal edges are repeated
  E = sort( ...
    [F(:,1) (:,2); ...
     F(:,2) (:,3); ...
     F(:,3) (:,1)]')';
  % remove duplicates
  E = unique(E,'rows');
  writeTGF(output_name,V,E);

end
