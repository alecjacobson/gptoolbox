if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 33a931d0191344e56007a0f709d6ce38af2e1ebd
)
FetchContent_MakeAvailable(libigl)
