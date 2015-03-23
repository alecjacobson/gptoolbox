MEXOPTS={'-v','-largeArrayDims','-DMEX'};
MSSE42='CXXFLAGS=\$CXXFLAGS -msse4.2';
EIGEN_INC='-I/opt/local/include/eigen3';
LIBIGL_INC=sprintf('-I%s/include',path_to_libigl);
LIBIGL_LIB=sprintf('-L%s/lib -ligl',path_to_libigl);
LIBIGL_LIBMATLAB='-liglembree';
LIBIGL_LIBEMBREE='-liglmatlab';
EMBREE=[path_to_libigl '/external/embree'];
EMBREE_INC=sprintf('-I%s -I%s/embree/',EMBREE,EMBREE);
EMBREE_LIB=sprintf('-L%s/build -lembree -lsys',EMBREE);

% bone_visible
mex(MEXOPTS{:},'-o','bone_visible_mex','bone_visible.cpp');
% ray_mesh_intersect
mex(MEXOPTS{:},'-o','ray_mesh_intersect','ray_mesh_intersect.cpp');
% bone_visible_embree
mex( ...
  MEXOPTS{:},...
  MSSE42, ...
  LIBIGL_INC,LIBIGL_LIB,LIBIGL_LIBEMBREE,LIBIGL_LIBMATLAB,EIGEN_INC, ...
  EMBREE_INC, EMBREE_LIB, ...
  '-o','bone_visible_embree', ...
  'bone_visible_embree.cpp');
