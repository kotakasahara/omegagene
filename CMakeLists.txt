CMAKE_MINIMUM_REQUIRED( VERSION 3.1.0 )

# CMake modules/macros are in a subdirectory to keep this file cleaner
# This needs to be set before PROJECT() in order to pick up toolchain files
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build )
PROJECT( CELESTE )
INCLUDE(SetCMakeMessageColors)

# SET(CMAKE_VERBOSE_MAKEFILE ON)
SET( _CMAKE_PROJECT_INITIALS "celeste" CACHE INTERNAL "Project Name Initials" )

MESSAGE( "\n${BoldCyan}    **** WELCOME TO ${CMAKE_PROJECT_NAME} (${_CMAKE_PROJECT_INITIALS})! ****${ColourReset}\n" )

INCLUDE(CheckCXXCompiler)
INCLUDE(CheckCUDA)
#INCLUDE(CheckOpenMP)

# INCLUDE_DIRECTORIES(include)

# Build project in the following subdirectories
ADD_SUBDIRECTORY(src)


# Remove all files not tracked by git
ADD_CUSTOM_TARGET(gitclean COMMAND git clean -d -f -x)
MESSAGE("\n${Magenta}    To clean out all files not tracked by git,\n    simply run 'make gitclean' from the\n    project ROOT directory.${ColourReset}\n")
