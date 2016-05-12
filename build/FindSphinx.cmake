find_program(SPHINX_EXECUTABLE NAMES sphinx-build sphinx-build2
    HINTS
    $ENV{SPHINX_DIR}
    PATH_SUFFIXES bin
    DOC "Sphinx documentation generator"
)
mark_as_advanced(SPHINX_EXECUTABLE)

# Detect Sphinx version
if (SPHINX_EXECUTABLE AND NOT DEFINED SPHINX_EXECUTABLE_VERSION)
    execute_process(
        COMMAND ${SPHINX_EXECUTABLE} --version
        OUTPUT_VARIABLE SPHINX_VERSION_OUTPUT_VARIABLE
        RESULT_VARIABLE SPHINX_VERSION_RESULT_VARIABLE
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    string(REGEX REPLACE "Sphinx \\([^)]*\\) ([^ ]+)" "\\1" SPHINX_EXECUTABLE_VERSION "${SPHINX_VERSION_OUTPUT_VARIABLE}")
    set(SPHINX_EXECUTABLE_VERSION "${SPHINX_EXECUTABLE_VERSION}" CACHE INTERNAL "Version of ${SPHINX_EXECUTABLE}")
endif()

set(_find_deps_options)
#if (Sphinx_FIND_QUIETLY)
#    set(_find_deps_options QUIET)
#endif()

include(FindPythonModule)
find_python_module(pygments ${_find_deps_options})
if (PYTHONMODULE_PYGMENTS)
    set(Sphinx_pygments_FOUND 1)
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Sphinx
    REQUIRED_VARS SPHINX_EXECUTABLE
    VERSION_VAR SPHINX_EXECUTABLE_VERSION
    HANDLE_COMPONENTS
)
