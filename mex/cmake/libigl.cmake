if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 1b07ebcf9e18c33ba342fd5c945dedbda055df78
)
FetchContent_MakeAvailable(libigl)
