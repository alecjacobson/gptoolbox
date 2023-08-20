cmake_minimum_required(VERSION 3.1)
set(CMAKE_BUILD_TYPE_INIT Release)
project(gptoolbox_mex)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# templated lambdas
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

# Move output to this dir (gptoolbox/mex)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

# https://stackoverflow.com/q/46724267/148668
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

# From libigl/cmake/libigl.cmake libigl/cmake/LibiglWindows.cmake
if(MSVC)
  if("${MSVC_RUNTIME}" STREQUAL "")
    set(MSVC_RUNTIME "static")
  endif()
  if(${MSVC_RUNTIME} STREQUAL "static")
    message(STATUS "MSVC -> forcing use of statically-linked runtime.")
    foreach(config ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER ${config} config)
      string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
      string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
    endforeach()
    string(REPLACE "/MDd" "/MTd" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
  else()
    message(STATUS "MSVC -> forcing use of dynamically-linked runtime.")
    foreach(config ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER ${config} config)
      string(REPLACE /MT /MD CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
      string(REPLACE /MT /MD CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
    endforeach()
    string(REPLACE "/MTd" "/MDd" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
  endif()
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(GPTOOLBOX_MEX_ROOT "${PROJECT_SOURCE_DIR}")
set(GPTOOLBOX_MEX_EXTERNAL "${GPTOOLBOX_MEX_ROOT}/external")
include(gptoolbox_download_project)

# libigl
option(LIBIGL_COPYLEFT_CGAL             "Use CGAL"           ON)
option(LIBIGL_EMBREE           "Use Embree"         ON)
# just to get it to download tetgen
option(LIBIGL_COPYLEFT_TETGEN           "Use Tegen"          ON)
option(LIBIGL_RESTRICTED_TRIANGLE         "Use Triangle"       ON)
option(LIBIGL_RESTRICTED_MATLAB           "Use Matlab"         ON)
option(LIBIGL_XML              "Use XML"            ON)
option(WITH_ELTOPO                  "Use El Topo"        ON)
option(MATLAB_PROXY                 "Download Matlab Proxy for compilation only" OFF)
if( UNIX AND NOT APPLE)
  # Static boost libs are causing problems on linux
  option(CGAL_Boost_USE_STATIC_LIBS   "Use Static  Boost libs" OFF)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif()

# Find matlab
if(MATLAB_PROXY)
  set(Matlab_ROOT_DIR "${GPTOOLBOX_MEX_ROOT}/external/matlab")
  gptoolbox_download_matlab()
endif()

find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY)

# Download and find libigl
include(libigl)

if(WITH_ELTOPO)
  # Download and find eltopo
  find_package(BLAS REQUIRED )
  gptoolbox_download_eltopo()
  set(LIBELTOPO_INCLUDE_DIR "${GPTOOLBOX_MEX_EXTERNAL}/eltopo/eltopo3d")
  add_subdirectory(${LIBELTOPO_INCLUDE_DIR} libeltopo)
endif()

set( CMAKE_SHARED_LIBRARY_SUFFIX ".${Matlab_MEX_EXTENSION}" )  # set suffix to .mexa64
set( CMAKE_SHARED_LIBRARY_PREFIX )      # remove the "lib" prefix
set( CMAKE_SHARED_LIBRARY_CXX_FLAGS )       # remove the -fPIC option. mex does not accept the "-fPIC" option
set( CMAKE_POSITION_INDEPENDENT_CODE ON)
set( CMAKE_INSTALL_RPATH "\$ORIGIN" )

function(compile_each I_MEX_SOURCES I_LIBS )
  foreach(MEX_SOURCE ${I_MEX_SOURCES})
    message(STATUS "MEX_SOURCE: ${MEX_SOURCE}")
    get_filename_component(MEX_OUTPUT ${MEX_SOURCE} NAME_WE)
    message(STATUS "MEX_OUTPUT: ${MEX_OUTPUT}")
    add_library(${MEX_OUTPUT} SHARED ${MEX_SOURCE})
    # Not sure if this should be (WIN32) or (MSVC OR MSYS OR MINGW)
    # https://stackoverflow.com/a/40152725/148668
    if(WIN32)
      # https://stackoverflow.com/a/11095194/148668
      set_target_properties(${MEX_OUTPUT} PROPERTIES LINK_FLAGS "/export:mexFunction")
    endif()
    target_link_libraries(${MEX_OUTPUT} ${I_LIBS})
    target_include_directories(${MEX_OUTPUT} PUBLIC ${Matlab_INCLUDE_DIRS})
    list(APPEND MEX_OUTPUTS ${MEX_OUTPUT})
  endforeach()
endfunction()

#file( GLOB MEX_SOURCES *.cpp )
list(APPEND CORE_LIBS ${Matlab_LIBRARIES})

#list(APPEND MEX_SOURCES "impaste.cpp")

# LIBIGL
add_definitions("-DMEX")
compile_each("\
aabb.cpp;\
angle_derivatives.cpp;\
bone_visible.cpp;\
blue_noise.cpp;\
collapse_small_triangles.cpp;\
decimate_libigl.cpp;\
dual_laplacian.cpp;\
exact_geodesic.cpp;\
fast_sparse.cpp;\
fit_rotations_mex.cpp;\
fit_cubic_bezier.cpp;\
icp.cpp;\
in_element_aabb.cpp;\
psd_project_rows.cpp;\
principal_curvature.cpp;\
readMSH.cpp;\
read_triangle_mesh.cpp;\
segment_graph.cpp;\
signed_distance.cpp;\
simplify_polyhedron.cpp;\
slim.cpp;\
solid_angle.cpp;\
"
  "${CORE_LIBS};igl::core;igl_restricted::matlab")

if(WITH_ELTOPO)
  compile_each(
    "eltopo.cpp" 
    "${CORE_LIBS};igl::core;igl_restricted::matlab;libeltopo;${BLAS_LIBRARIES}"
    )
endif()

if(LIBIGL_COPYLEFT_TETGEN)
  compile_each("\
tetrahedralize.cpp;\
"
  "${CORE_LIBS};igl::core;igl_restricted::matlab;igl_copyleft::tetgen")
endif()

if(LIBIGL_RESTRICTED_TRIANGLE)
  if(LIBIGL_COPYLEFT_CGAL)
    add_definitions("-DWITH_CGAL")
    compile_each("\
triangulate.cpp;\
"
    "${CORE_LIBS};igl::core;igl_restricted::matlab;igl_restricted::triangle;igl_copyleft::cgal")
  else()
    compile_each("\
triangulate.cpp;\
"
    "${CORE_LIBS};igl::core;igl_restricted::matlab;igl_restricted::triangle")
  endif()
endif()

if(LIBIGL_EMBREE)
  add_definitions("-DWITH_EMBREE")
  compile_each("\
ambient_occlusion.cpp;\
bone_visible_embree.cpp;\
ray_mesh_intersect.cpp;\
reorient_facets.cpp;\
"
    "${CORE_LIBS};igl::core;igl_restricted::matlab;igl::embree")
endif()

if(LIBIGL_EMBREE)
  if(LIBIGL_COPYLEFT_CGAL)
    compile_each("winding_number.cpp" 
      "${CORE_LIBS};igl::core;igl_restricted::matlab;igl::embree;igl_copyleft::cgal")
  else()
    compile_each("winding_number.cpp" 
      "${CORE_LIBS};igl::core;igl_restricted::matlab;igl::embree")
  endif()
else()
  compile_each("winding_number.cpp" 
    "${CORE_LIBS};igl::core;igl_restricted::matlab")
endif()



if(LIBIGL_COPYLEFT_CGAL)
  add_definitions("-DWITH_CGAL")
#segment_segment_intersection.cpp;\
  compile_each("\
box_intersect.cpp;\
decimate_cgal.cpp;\
form_factor.cpp;\
intersect_other.cpp;\
mesh_boolean.cpp;\
outer_hull.cpp;\
point_mesh_squared_distance.cpp;\
selfintersect.cpp;\
signed_distance_isosurface.cpp;\
snap_rounding.cpp;\
trim_with_solid.cpp;\
upper_envelope.cpp;\
wire_mesh.cpp;\
"
    "${CORE_LIBS};igl::core;igl_restricted::matlab;igl_copyleft::cgal;${CORK_LIBS}")
  if(LIBIGL_XML)
    compile_each(
      "read_mesh_from_xml.cpp"
      "${CORE_LIBS};igl::core;igl_restricted::matlab;igl_copyleft::cgal;igl::xml")
  endif()
endif()

if(APPLE)
  find_library(FOUNDATION_LIBRARY Foundation)
  find_library(APPKIT_LIBRARY AppKit)
  mark_as_advanced (FOUNDATION_LIBRARY APPKIT_LIBRARY)
  message(STATUS "MEX_SOURCE: impaste.cpp")
  set(MEX_OUTPUT "impaste")
  add_library(${MEX_OUTPUT} SHARED "impaste.cpp" "paste.mm")
  target_link_libraries(${MEX_OUTPUT} 
    ${CORE_LIBS} ${FOUNDATION_LIBRARY} ${APPKIT_LIBRARY})
  target_include_directories(${MEX_OUTPUT} PUBLIC ${Matlab_INCLUDE_DIRS})
endif(APPLE)
