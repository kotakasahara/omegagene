CMAKE_MINIMUM_REQUIRED( VERSION 3.1.0 )

# CMake modules/macros are in a subdirectory to keep this file cleaner
# This needs to be set before PROJECT() in order to pick up toolchain files
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build )
PROJECT( CELESTE )

# if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR} )
#     MESSAGE( FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt." )
# endif()

# SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

INCLUDE(SetCMakeMessageColors)

# SET(CMAKE_VERBOSE_MAKEFILE ON)
SET( _CMAKE_PROJECT_INITIALS "celeste" CACHE INTERNAL "Project Name Initials" )

MESSAGE( "\n${BoldCyan}    **** WELCOME TO ${CMAKE_PROJECT_NAME} (${_CMAKE_PROJECT_INITIALS})! ****${ColourReset}\n" )

INCLUDE(CheckCXXCompiler)
INCLUDE(CheckCelesteBuildOptions)

IF(CELESTE_GPU OR CELESTE_GPUECP)
	INCLUDE(CheckCUDA)
ENDIF()
#INCLUDE(CheckOpenMP)

# INCLUDE_DIRECTORIES(include)

# Build project in the following subdirectories
ADD_SUBDIRECTORY(src)


# Remove all files not tracked by git
ADD_CUSTOM_TARGET(gitclean COMMAND git clean -d -f -x)
MESSAGE("\n${Magenta}    To clean out all files not tracked by git,\n    simply run 'make gitclean' from the\n    project ROOT directory.${ColourReset}\n")
