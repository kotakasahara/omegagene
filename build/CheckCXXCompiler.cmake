# Use C++11
INCLUDE(CheckCXXCompilerFlag)
SET( MINIMUM_GCC_VERSION_REQUIRED 4.8 )
SET( MINIMUM_CLANG_VERSION_REQUIRED 3.1 )
SET( MINIMUM_ICC_VERSION_REQUIRED 15.0 )

MACRO(CHECK_COMPILER_VERSION compiler_name minimum_version)
    MESSAGE( "-- Checking for C++ compiler version" )
    IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${minimum_version})
        MESSAGE( FATAL_ERROR "${BoldRed}${compiler_name} version must be >= ${MINIMUM_GCC_VERSION_REQUIRED}!${ColourReset}" )
    ENDIF()
    MESSAGE( "-- Checking for C++ compiler version -> ${CMAKE_CXX_COMPILER_VERSION} (Success)" )
ENDMACRO(CHECK_COMPILER_VERSION)


CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
IF (COMPILER_SUPPORTS_CXX14)
    MESSAGE(STATUS "The installed ${CMAKE_CXX_COMPILER_ID} compiler has C++14 support!")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
ELSEIF (COMPILER_SUPPORTS_CXX11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSE()
    MESSAGE(FATAL_ERROR "${BoldRed}This project requires C++11 and the installed ${CMAKE_CXX_COMPILER_ID} compiler has no C++11 support. Please install a modern version of ${CMAKE_CXX_COMPILER_ID}!${ColourReset}")
ENDIF()

# On Linux machines, clang uses gcc's libstdc++ library by default
# External libraries built using libstdc++ will need to be rebuilt with libc++ though - see http://stackoverflow.com/questions/8454329/why-cant-clang-with-libc-in-c0x-mode-link-this-boostprogram-options-examp
CHECK_CXX_COMPILER_FLAG("-stdlib=libc++" COMPILER_SUPPORTS_LIBCXX)
IF (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND COMPILER_SUPPORTS_LIBCXX)
    MESSAGE( STATUS "libc++ found; explicitly compiling with -stdlib=libc++..." )
    MESSAGE( WARNING "${UnderlinedRed}External dependencies may need to be rebuilt using libc++ to avoid linker errors during project compilation!${ColourReset}" )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
ENDIF()

# Older versions of compilers may support C++11, but not include C++11 standards-compliant libraries
IF (CMAKE_COMPILER_IS_GNUCC)
    CHECK_COMPILER_VERSION(GCC MINIMUM_GCC_VERSION_REQUIRED)
ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    CHECK_COMPILER_VERSION(ICC MINIMUM_ICC_VERSION_REQUIRED)
ELSEIF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    CHECK_COMPILER_VERSION(Clang MINIMUM_CLANG_VERSION_REQUIRED)
ENDIF()

# Enable compiler warnings
IF (MSVC)
    IF (CMAKE_CXX_FLAGS MATCHES "/W[0-4]") # Force to always compile with W4
        STRING(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    ELSE()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    ENDIF()
ELSEIF (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
ENDIF()
