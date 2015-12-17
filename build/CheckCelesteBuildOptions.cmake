MACRO(SET_BUILD_FLAGS flag_name description_string default_value)
    SET(${flag_name} ${default_value})  # Remove the flag from CMakeCache.txt so that CMake does not cache flag values
    OPTION(${flag_name} ${description_string} ${default_value})
    MESSAGE("[BUILD OPTIONS] ${flag_name}=${${flag_name}}")
ENDMACRO(SET_BUILD_FLAGS)

SET_BUILD_FLAGS(CELESTE_WO_NS    "Without neighbor-search\; all pairs of atoms are evaluated"                    0)
SET_BUILD_FLAGS(CELESTE_GPU      "Neighbor-search and pairwise calculation are done on GPU"                      0)
SET_BUILD_FLAGS(CELESTE_GPUECP   "Neighbor-search is done on CPU, pairwise potential calculation is done on GPU" 0)

IF ((CELESTE_WO_NS AND CELESTE_GPU) OR (CELESTE_GPU AND CELESTE_GPUECP) OR (CELESTE_GPUECP AND CELESTE_WO_NS))
    MESSAGE( FATAL_ERROR "${BoldRed}Too many build configs chosen; please choose only one of [ CELESTE_WO_NS, CELESTE_GPU, CELESTE_GPUECP ]${ColourReset}")
ENDIF()
