if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG db07a47bec950081309fdca775fcfd86f2f04af2
)
FetchContent_MakeAvailable(libigl)
