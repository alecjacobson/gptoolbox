################################################################################
include(DownloadProject)

# Cribbed from LibiglDownloadExternal.cmake

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
	set(GPTOOLBOX_DOWNLOAD_PROJECT_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
	set(GPTOOLBOX_DOWNLOAD_PROJECT_OPTIONS "")
endif()

# Shortcut function
function(gptoolbox_download_project name)
	download_project(
		PROJ         ${name}
                SOURCE_DIR   ${GPTOOLBOX_MEX_EXTERNAL}/${name}
		DOWNLOAD_DIR ${GPTOOLBOX_MEX_EXTERNAL}/.cache/${name}
		QUIET
		${GPTOOLBOX_DOWNLOAD_PROJECT_OPTIONS}
		${ARGN}
	)
endfunction()

################################################################################

## Matlab Proxy
function(gptoolbox_download_matlab)
	gptoolbox_download_project(matlab
		GIT_REPOSITORY https://github.com/alecjacobson/matlab.git
		GIT_TAG        883d417a99fcb8ead89387cee243e51a92864019
	)
endfunction()

## LIBIGL
function(gptoolbox_download_libigl)
	gptoolbox_download_project(libigl
		GIT_REPOSITORY https://github.com/libigl/libigl.git
		GIT_TAG        a5a501f3a7c010ca5f3c5da81e471f4d33cf57f0
	)
endfunction()

## EL TOPO
function(gptoolbox_download_eltopo)
	gptoolbox_download_project(eltopo
		GIT_REPOSITORY https://github.com/alecjacobson/eltopo.git
		GIT_TAG        e4c2c9a3f01b89bcafb3afa203d2bad887f4133b
	)
endfunction()
