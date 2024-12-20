if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 667101084ac342f1f700299349452143069cbb4a
)
FetchContent_MakeAvailable(libigl)
