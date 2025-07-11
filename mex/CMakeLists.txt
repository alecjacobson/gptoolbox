cmake_minimum_required(VERSION 3.16.0)
set(CMAKE_BUILD_TYPE_INIT Release)
project(gptoolbox_mex)
cmake_policy(SET CMP0074 NEW)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# templated lambdas
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

# Move output to this dir (gptoolbox/mex)
if(MSVC)
  # Solution for multitarget-generators like VS2019 From: https://stackoverflow.com/questions/47175912/using-cmake-how-to-stop-the-debug-and-release-subdirectories
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $<1:${PROJECT_SOURCE_DIR}>)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY $<1:${PROJECT_SOURCE_DIR}>)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY $<1:${PROJECT_SOURCE_DIR}>)
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
endif()

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

# libigl
option(LIBIGL_COPYLEFT_CGAL       "Use CGAL"           ON)
option(LIBIGL_EMBREE              "Use Embree"         ON)
option(LIBIGL_PREDICATES          "Use Predicates"         ON)
option(LIBIGL_COPYLEFT_TETGEN     "Use Tetgen"          ON)
option(LIBIGL_RESTRICTED_TRIANGLE "Use Triangle"       ON)
option(LIBIGL_RESTRICTED_MATLAB   "Use Matlab"         ON)
option(LIBIGL_XML                 "Use XML"            ON)
option(WITH_ELTOPO                "Use El Topo"        ON)
option(MATLAB_PROXY               "Download Matlab Proxy for compilation only" OFF)
if( UNIX AND NOT APPLE)
  # Static boost libs are causing problems on linux
  option(CGAL_Boost_USE_STATIC_LIBS   "Use Static  Boost libs" OFF)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif()

# git
find_package(Git REQUIRED)

# ---[ setup vcpkg toolchain ]---
# vcpkg defines a set of variables if toolchain setup was successfull
if(NOT DEFINED Z_VCPKG_EXECUTABLE)
  set(Z_VCPKG_ROOT_DIR "${GPTOOLBOX_MEX_EXTERNAL}/vcpkg")
  message(STATUS "vcpkg was not found in system and will be installed to ${Z_VCPKG_ROOT_DIR} ...")
  if(NOT EXISTS ${Z_VCPKG_ROOT_DIR})
    execute_process(
      COMMAND             "${GIT_EXECUTABLE}" clone https://github.com/microsoft/vcpkg "${Z_VCPKG_ROOT_DIR}"
      COMMAND_ERROR_IS_FATAL ANY
    )
  endif()
  
  # CMAKE_EXECUTABLE_SUFFIX is not yet defined
  if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
    set(Z_VCPKG_EXECUTABLE "${Z_VCPKG_ROOT_DIR}/vcpkg.exe")
    set(Z_VCPKG_BOOTSTRAP_SCRIPT "${Z_VCPKG_ROOT_DIR}/bootstrap-vcpkg.bat")
  else()
    set(Z_VCPKG_EXECUTABLE "${Z_VCPKG_ROOT_DIR}/vcpkg")
    set(Z_VCPKG_BOOTSTRAP_SCRIPT "${Z_VCPKG_ROOT_DIR}/bootstrap-vcpkg.sh")
  endif()


  # run bootstrap once
  if(NOT EXISTS ${Z_VCPKG_EXECUTABLE})
    execute_process(
      COMMAND ${Z_VCPKG_BOOTSTRAP_SCRIPT} -disableMetrics
      COMMAND_ERROR_IS_FATAL ANY
    )
  endif()

  set(VCPKG_ARGS "--clean-after-build")
endif() #vcpkg
# toolchain definition
set(CMAKE_TOOLCHAIN_FILE "${Z_VCPKG_ROOT_DIR}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
include("${Z_VCPKG_ROOT_DIR}/scripts/buildsystems/vcpkg.cmake")

# Find matlab
if(MATLAB_PROXY)
  FetchContent_Declare(
    Matlab
    GIT_REPOSITORY https://github.com/alecjacobson/matlab.git
    GIT_TAG        883d417a99fcb8ead89387cee243e51a92864019
  )
  FetchContent_GetProperties(Matlab)
  if(NOT Matlab_POPULATED)
    FetchContent_Populate(Matlab)
    set(Matlab_ROOT_DIR "${Matlab_SOURCE_DIR}")
  endif()
endif()
find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY)

message(STATUS "Downloading and finding libigl...")
# Download and find libigl
# NOTE: should be replaced with vcpkg version of libigl for more consistency as soon apdated to 2.4.0 (and https://github.com/microsoft/vcpkg/pull/26029 colsed) there
#       find_package(LIBIGL REQUIRED)
include(libigl)


if(WITH_ELTOPO)
  # Download and find eltopo
  find_package(BLAS REQUIRED )
  FetchContent_Declare(
    eltopo
    GIT_REPOSITORY https://github.com/alecjacobson/eltopo.git
    GIT_TAG        3cf988aae467392d1066bdde0be3f54487a0ba2c
  )
  FetchContent_GetProperties(eltopo)
  if(NOT eltopo_POPULATED)
    FetchContent_Populate(eltopo)
    # Point to the directory containing the CMakeLists.txt of libeltopo
    add_subdirectory(${eltopo_SOURCE_DIR}/eltopo3d ${eltopo_BINARY_DIR})
  endif()
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
isolines.cpp;\
in_element_aabb.cpp;\
psd_project_rows.cpp;\
principal_curvature.cpp;\
readMSH.cpp;\
read_triangle_mesh.cpp;\
segment_graph.cpp;\
signed_distance.cpp;\
simplify_polyhedron.cpp;\
slim.cpp;\
split_nonmanifold.cpp;\
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
  compile_each("\
refine_triangulation.cpp;\
"
  "${CORE_LIBS};igl::core;igl_restricted::matlab;igl_restricted::triangle")
endif()

if(LIBIGL_PREDICATES)
  add_definitions("-DWITH_PREDICATES")
  compile_each("\
orient3d.cpp;\
"
    "${CORE_LIBS};igl::core;igl_restricted::matlab;igl::predicates")
endif()

if(LIBIGL_EMBREE)
  add_definitions("-DWITH_EMBREE")
  compile_each("\
ambient_occlusion.cpp;\
bone_visible_embree.cpp;\
ray_mesh_intersect.cpp;\
ray_mesh_intersect_all.cpp;\
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
  compile_each("\
box_intersect.cpp;\
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
