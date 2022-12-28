if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 3ea7f9480967fcf6bf02ce9b993c0ea6d2fc45f6
)
FetchContent_MakeAvailable(libigl)
