function [V,T,F] = tetgen(SV,SF,IV,allow_resampling)
  % TETGEN
  % [V,T,F] = tetgen(SV,SF,IV)
  %
  % Call tetgen to construct a tetrahedral volume mesh with in a given triangle
  % mesh with optional internal contrained vertices.
  %
  % Inputs:
  %   SV  list of surface vertex positions of exterior mesh, # vertices by 3
  %   SF  list of surface face indices of exterior triangle mesh, # faces by 3
  %   IV  list of internal vertex positions, # internal vertice by 3
  %   allow_resampling  allow resampling on the surface given [false]
  % Outputs:
  %   V  list of tetrahedra vertices
  %   T  list of tetrahedra indices
  %   F  list of faces of 3D volume mesh
  %

  % determine if internal constraint vertices are present
  internal_constraints = false;
  if(exist('IV','var') && prod(size(IV))>0)
    internal_constraints = true;
    assert(max(size(setdiff(SV,IV,'rows')) == size(SV)))
  end

  % get a temporary file name prefix
  prefix = tempname;
  off_filename = [prefix '.off'];
  writeOFF(off_filename,SV,SF);

  % if there are internal constraint vertices then print them to a .node file
  if(internal_constraints)
    inode_filename = [prefix '.a.node'];
    writeNODE(inode_filename,IV);
  end

  path_to_tetgen = '/usr/local/bin/tetgen';
  % graded: -q100, very-fine:-q1
  flags = '-Cp -q100 ';
  if(internal_constraints)
    flags = [flags ' -i'];
  end
  %if(~exist('allow_resampling') || ~allow_resampling)
  %  flags = [flags ' -Y' '-V'];
  %end
  % call tetgen
  command = [path_to_tetgen ' ' flags ' ' off_filename];
  fprintf(command);
  [status, result] = system(command);
  status
  result

  % tetgen always writes output to file:
  %   xxxx.1.ele  tetrahedra
  %   xxxx.1.node tetrahedra vertices
  %   xxxx.1.face  surface faces
  ele_filename = [prefix '.1.ele'];
  face_filename = [prefix '.1.face'];
  node_filename = [prefix '.1.node'];

  F = readFACE(face_filename);
  % reverse faces because tetgen uses backwards order
  F = fliplr(F);
  T = readELE(ele_filename);
  V = readNODE(node_filename);


  delete(off_filename);
  if(internal_constraints)
    delete(inode_filename);
  end
  delete(ele_filename);
  delete(face_filename);
  delete(node_filename);

end
