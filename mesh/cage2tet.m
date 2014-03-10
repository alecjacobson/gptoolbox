function cage2tet( ...
  cage_name, ...
  surface_name, ...
  skeleton_name, ...
  samples_per_edge, ...
  output_name, ...
  allow_resampling)
  % CAGE2TET Script to make a 3d tet mesh from a 2d surface cage and a 2d
  % surface residing in the cage. Faces of the resulting .mesh file are the
  % faces of the original internal surface though they may not correspond to
  % faces of tets
  %
  % cage2tet( ...
  %   cage_name, surface_name, skeleton_name, ...
  %   samples_per_edge, output_name, allow_resampling)
  %
  % Inputs:
  %   cage_name  path to the file containing cage vertices and faces, should be
  %     off or obj
  %   surface_name  path to the file containing the internal surface vertices
  %     and faces
  %   skeleton_name  path to .tgf file containing skeleton
  %   samples_per_edge  number of samples per skeleton edge
  %   output_name  path to file to be written with tet mesh
  %   allow_resampling  whether to allow resampling on cage surface {false}
  %

  if(~exist('allow_resampling','var'))
    allow_resampling = false;
  end

  % read input mesh
  [cageV, cageF] = load_mesh(cage_name);

  % read input mesh
  [surfV,surfF] = load_mesh(surface_name);

  % read skelton
  if(isempty(skeleton_name))
    skelV = [];
  else
    [skelV,skelF] = readTGF(skeleton_name);
    skelV = sample_edges(skelV,skelF,samples_per_edge);
  end

  assert(max(size(setdiff(cageV,surfV,'rows')) == size(cageV)))
  assert(all(size(unique(surfV,'rows')) == size(surfV)))
  assert(all(size(unique(cageV,'rows')) == size(cageV)))

  [V,T,F] = tetgen([cageV;surfV],cageF,[skelV],allow_resampling);
  offset = find( ...
    V(:,1)==surfV(1,1) & ...
    V(:,2)==surfV(1,2) & ...
    V(:,2)==surfV(1,2) ...
    ,1)-1;
  indices = 1:size(V,1);
  indices(1:size(cageV,1)) = size(surfV,1) + (1:size(cageV,1));
  indices(size(cageV,1) + (1:size(surfV,1))) = 1:size(surfV,1);
  mT = indices(T);
  mV = zeros(size(V));
  mV(indices,:) = V;
  mF = surfF;

  %% write tet mesh to file
  writeMESH(output_name,mV,mT,mF);

  % subplot(1,2,1);tsurf(mF,mV);
  % subplot(1,2,2);tetramesh(mT,mV);
end
