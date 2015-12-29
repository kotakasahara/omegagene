FIND_PACKAGE(Git QUIET)
IF (NOT GIT_FOUND)
    MESSAGE( FATAL_ERROR "${BoldRed}Did not find git installed!${ColourReset}" )
ENDIF()

IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    MESSAGE( FATAL_ERROR "${BoldRed}/.git/ not existent; cannot determine CELESTE_BUILD_VERSION!${ColourReset}" )
ENDIF()

# Get project version from git tag
EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} describe --long --tags --dirty --always
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE CELESTE_BUILD_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# If the git tree is dirty, use current timestamp; else use the latest commit's timestamp
IF (CELESTE_BUILD_VERSION MATCHES "dirty$")
    STRING(TIMESTAMP CELESTE_BUILD_TIMESTAMP "%Y-%m-%d %H:%M:%S")
ELSE()
    EXECUTE_PROCESS(
        COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:"%cd" --date=iso
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE CELESTE_BUILD_TIMESTAMP
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    STRING(REGEX REPLACE "\"" "" CELESTE_BUILD_TIMESTAMP ${CELESTE_BUILD_TIMESTAMP}) # Remove quotes
ENDIF()

MESSAGE( "${BoldBlue}    Build Version   : ${CELESTE_BUILD_VERSION}\n    Build Timestamp : ${CELESTE_BUILD_TIMESTAMP}${ColourReset}\n" )
ADD_DEFINITIONS( -DMACRO_BUILD_VERSION="${CELESTE_BUILD_VERSION}" -DMACRO_BUILD_TIMESTAMP="${CELESTE_BUILD_TIMESTAMP}" )
