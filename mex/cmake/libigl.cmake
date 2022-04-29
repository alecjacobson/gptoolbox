if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 45a4f79b5ec72327c59758fc87f6fb72a9edffc8
)
FetchContent_MakeAvailable(libigl)
