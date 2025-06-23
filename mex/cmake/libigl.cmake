if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 0e360d5250b5c24b5c483d850c42cd92e0a09f53
)
FetchContent_MakeAvailable(libigl)
