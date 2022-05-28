if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 30019acb17197abb6f3fe4c310909c817f2e85cc
)
FetchContent_MakeAvailable(libigl)
