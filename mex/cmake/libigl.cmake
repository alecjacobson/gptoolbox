if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 73a2e0de9b0bcfac83dc5bb661049c516a73650f
)
FetchContent_MakeAvailable(libigl)
