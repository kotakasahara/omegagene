CMAKE_MINIMUM_REQUIRED( VERSION 3.1.0 )

LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../build )
PROJECT( "CELESTE DOCS" )

include(SetCMakeMessageColors)
include(FindSphinx)

find_package(Sphinx REQUIRED)

if(NOT DEFINED SPHINX_THEME)
    set(SPHINX_THEME default)
endif()

if(NOT DEFINED SPHINX_THEME_DIR)
    set(SPHINX_THEME_DIR)
endif()

# source directory
set(SPHINX_SOURCES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")

# configured documentation tools and intermediate build results
set(SPHINX_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${SPHINX_BUILD_DIR}/doctrees")


add_custom_target(html ALL
	${SPHINX_EXECUTABLE} -q -b html -c "${SPHINX_SOURCES_DIR}" -d "${SPHINX_CACHE_DIR}" "${SPHINX_SOURCES_DIR}" "${SPHINX_BUILD_DIR}/html"
    COMMENT "Building HTML documentation with Sphinx")

add_custom_target(latexpdf ALL
	${SPHINX_EXECUTABLE} -q -b latexpdf -c "${SPHINX_SOURCES_DIR}" -d "${SPHINX_CACHE_DIR}" "${SPHINX_SOURCES_DIR}" "${SPHINX_BUILD_DIR}/latexpdf"
    COMMENT "Building PDF documentation with Sphinx")

add_custom_target(latexpdfja ALL
	${SPHINX_EXECUTABLE} -q -b latexpdfja -c "${SPHINX_SOURCES_DIR}" -d "${SPHINX_CACHE_DIR}" "${SPHINX_SOURCES_DIR}" "${SPHINX_BUILD_DIR}/latexpdfja"
    COMMENT "Building PDF documentation with Sphinx")

add_custom_target(manpage ALL
    ${SPHINX_EXECUTABLE} -q -b man -c "${SPHINX_SOURCES_DIR}" -d "${SPHINX_CACHE_DIR}" "${SPHINX_SOURCES_DIR}" "${SPHINX_BUILD_DIR}/manpage"
    COMMENT "Building manpage documentation with Sphinx")

message("\n${Cyan}  After running 'make', please open in your browser:\n  ${SPHINX_BUILD_DIR}/html/index.html   ${ColourReset}\n")
