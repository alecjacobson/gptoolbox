if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG fe7e8281a4d92d422d8e516e155c3dd7cd26ab19
)
FetchContent_MakeAvailable(libigl)
