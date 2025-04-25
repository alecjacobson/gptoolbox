if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG 89267b4a80b1904de3f6f2812a2053e5e9332b7e
)
FetchContent_MakeAvailable(libigl)
