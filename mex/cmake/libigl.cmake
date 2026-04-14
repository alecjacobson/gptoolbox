if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 6000ccb70fdeb78376dcb5d2531e57a15d884aa0
)
FetchContent_MakeAvailable(libigl)
