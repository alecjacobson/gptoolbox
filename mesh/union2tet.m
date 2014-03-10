function union2tet(surface_name,cage_name,union_name,tet_name)
  % UNION2TET Script to take 3 meshes, a surface, a partial cage that intersects the
  % surface and their volume-union's surface (probably from Blender), and use
  % tetgen to generate a tet mesh that has tets in the volume-union with faces
  % of the original surface .
  %
  % union2tet(surface_name,cage_name,union_name,tet_name)
  %
  % Inputs:
  %   surface_name  file containing original surface
  %   cage_name  file containing (partial) cage
  %   union_name  file containing union surface
  % Output:
  %   tet_name  file to be written with .mesh tet mesh
  %
  [SV,SF] = load_mesh(surface_name);
  [CV,CF] = load_mesh(cage_name);
  [UV,UF] = load_mesh(union_name);
  % combine union and surface so we have a list of faces for the union and the
  % original surface that look into the same vertex list
  [V,F] = combine(SV,SF,UV,UF);
  [T,V] = tetgen(V,F,setdiff(CV,V,'rows'),true);
  writeMESH(tet_name,V,T,SF);

  skeleton_name = strrep([tet_name '.tgf'],'.mesh','');
  % Find all edges in mesh, note internal edges are repeated
  CE = sort( ...
    [CF(:,1) CF(:,2); ...
    CF(:,2) CF(:,3); ...
    CF(:,3) CF(:,1)]')';
  % remove duplicates
  CE = unique(CE,'rows');
  writeTGF(skeleton_name,CV,CE);
  

end
