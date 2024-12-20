if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 68f54c43ccebe3165c5cad5aa15003db72410989
)
FetchContent_MakeAvailable(libigl)
