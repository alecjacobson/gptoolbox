if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 9162fb7d797c2a55ce15e467dcd36e6d9e8ab245
)
FetchContent_MakeAvailable(libigl)
