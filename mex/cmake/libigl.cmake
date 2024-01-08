if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 7f7f0fe00709407bbaf7905d0827d89e2d750781
)
FetchContent_MakeAvailable(libigl)
