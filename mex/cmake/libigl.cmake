if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 7313592e16d8d3db393166d0d636d6106ffa890d
)
FetchContent_MakeAvailable(libigl)
