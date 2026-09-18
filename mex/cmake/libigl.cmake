if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG c678e8658b9c09954285bdf9687bd4b195ab1c97
)
FetchContent_MakeAvailable(libigl)
