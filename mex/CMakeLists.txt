cmake_minimum_required(VERSION 3.1)
set(CMAKE_BUILD_TYPE_INIT Release)
project(gptoolbox_mex)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Move output to this dir (gptoolbox/mex)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})


list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(GPTOOLBOX_MEX_ROOT "${PROJECT_SOURCE_DIR}")
set(GPTOOLBOX_MEX_EXTERNAL "${GPTOOLBOX_MEX_ROOT}/external")
include(gptoolbox_download_project)

# libigl
option(LIBIGL_WITH_CGAL             "Use CGAL"           ON)
option(LIBIGL_WITH_EMBREE           "Use Embree"         ON)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         ON)
option(WITH_ELTOPO                  "Use El Topo"        ON)

# Find matlab
find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY)
include_directories(${Matlab_INCLUDE_DIRS})

# Download and find libigl
gptoolbox_download_libigl()
find_package(LIBIGL REQUIRED)

if(WITH_ELTOPO)
  # Download and find eltopo
  find_package(BLAS REQUIRED )
  gptoolbox_download_eltopo()
  find_package(ELTOPO REQUIRED)
endif()

set( CMAKE_SHARED_LIBRARY_SUFFIX ".${Matlab_MEX_EXTENSION}" )  # set suffix to .mexa64
set( CMAKE_SHARED_LIBRARY_PREFIX )      # remove the "lib" prefix
set( CMAKE_SHARED_LIBRARY_CXX_FLAGS )       # remove the -fPIC option. mex does not accept the "-fPIC" option
set( CMAKE_POSITION_INDEPENDENT_CODE ON)
set( CMAKE_INSTALL_RPATH "\$ORIGIN" )

#file( GLOB MEX_SOURCES *.cpp )
set(MEX_SOURCES)
set(LIBS)
list(APPEND LIBS ${Matlab_LIBRARIES})

#list(APPEND MEX_SOURCES "impaste.cpp")
list(APPEND MEX_SOURCES "segment_graph.cpp")

# LIBIGL
add_definitions("-DMEX")
list(APPEND LIBS "igl::core" "igl::matlab")
list(APPEND MEX_SOURCES
  "bone_visible.cpp"
  "collapse_small_triangles.cpp"
  "decimate_libigl.cpp"
  "fit_rotations_mex.cpp"
  "in_element_aabb.cpp"
  "aabb.cpp"
  "read_triangle_mesh.cpp"
  "readMSH.cpp"
  "signed_distance.cpp"
  "simplify_polyhedron.cpp"
  "slim.cpp"
  "solid_angle.cpp"
  "winding_number.cpp")

if(WITH_ELTOPO)
  list(APPEND MEX_SOURCES 
    "eltopo.cpp"
    )
  include_directories(${LIBELTOPO_INCLUDE_DIR})
  list(APPEND LIBS ${LIBELTOPO_LIBRARY})
  list(APPEND LIBS ${BLAS_LIBRARIES})
endif()

if(LIBIGL_WITH_EMBREE)
  list(APPEND MEX_SOURCES 
    "ambient_occlusion.cpp"
    "bone_visible_embree.cpp"
    "ray_mesh_intersect.cpp"
    "reorient_facets.cpp")
  list(APPEND LIBS igl::embree)
  # needed winding number
  add_definitions("-DWITH_EMBREE")
endif()

if(LIBIGL_WITH_CORK)
  list(APPEND LIBS igl::cork)
  add_definitions("-DWITH_CORK")
endif()

if(LIBIGL_WITH_CGAL)
  list(APPEND MEX_SOURCES 
    "decimate_cgal.cpp"
    "intersect_other.cpp"
    "mesh_boolean.cpp"
    "outer_hull.cpp"
    "point_mesh_squared_distance.cpp"
    "selfintersect.cpp"
    "signed_distance_isosurface.cpp"
    "snap_rounding.cpp"
    "trim_with_solid.cpp"
    "wire_mesh.cpp"
    )
  list(APPEND LIBS igl::cgal)
  # needed winding number
  add_definitions("-DWITH_CGAL")
endif()

if(LIBIGL_WITH_CGAL AND LIBIGL_WITH_XML)
  list(APPEND MEX_SOURCES
    "read_mesh_from_xml.cpp"
    )
  list(APPEND LIBS igl::xml)
endif()


message(STATUS "Compiling: ${MEX_SOURCES}")

set(MEX_OUTPUTS)
foreach(MEX_SOURCE ${MEX_SOURCES})
  get_filename_component(MEX_OUTPUT ${MEX_SOURCE} NAME_WE)
  add_library(${MEX_OUTPUT} SHARED ${MEX_SOURCE})
  target_link_libraries(${MEX_OUTPUT} ${LIBS})
  list(APPEND MEX_OUTPUTS ${MEX_OUTPUT})
endforeach()


if(APPLE)
  find_library(FOUNDATION_LIBRARY Foundation)
  find_library(APPKIT_LIBRARY AppKit)
  mark_as_advanced (FOUNDATION_LIBRARY APPKIT_LIBRARY)
  set(MEX_OUTPUT "impaste")
  add_library(${MEX_OUTPUT} SHARED "impaste.cpp" "paste.mm")
  target_link_libraries(${MEX_OUTPUT} 
    ${LIBS} ${FOUNDATION_LIBRARY} ${APPKIT_LIBRARY})
ENDIF (APPLE)
