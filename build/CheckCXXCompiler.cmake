# Use C++11
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if (COMPILER_SUPPORTS_CXX14)
    MESSAGE(STATUS "The installed ${CMAKE_CXX_COMPILER_ID} has C++14 support!")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
elseif (COMPILER_SUPPORTS_CXX11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    MESSAGE(FATAL_ERROR "${BoldRed}This project requires C++11 and the installed ${CMAKE_CXX_COMPILER_ID} has no C++11 support. Please install a modern version of ${CMAKE_CXX_COMPILER_ID}!${ColourReset}")
endif()

# On Linux machines, clang uses gcc's libstdc++ library by default
# External libraries built using libstdc++ will need to be rebuilt with libc++ though - see http://stackoverflow.com/questions/8454329/why-cant-clang-with-libc-in-c0x-mode-link-this-boostprogram-options-examp
CHECK_CXX_COMPILER_FLAG("-stdlib=libc++" COMPILER_SUPPORTS_LIBCXX)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND COMPILER_SUPPORTS_LIBCXX)
    MESSAGE( STATUS "libc++ found; explicitly compiling with -stdlib=libc++..." )
    MESSAGE( WARNING "${UnderlinedRed}External dependencies may need to be rebuilt using libc++ to avoid linker errors during project compilation!${ColourReset}" )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

# Older versions of compilers may support C++11, but not include C++11 standards-compliant libraries
MESSAGE( "-- Checking for C++ Standard Library conformance" )
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
    MESSAGE( FATAL_ERROR "${BoldRed}GCC version must be >= 5.0; only GCC 5.0+ comes with a fully C++11-compliant standard library!${ColourReset}" )
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.1)
    MESSAGE( FATAL_ERROR "${BoldRed}Clang version must be >= 3.1; only Clang 3.1+ comes with a fully C++11-compliant standard library!${ColourReset}" )
endif()
MESSAGE( "-- Checking for C++ Standard Library conformance - Success" )

# Enable compiler warnings
if (MSVC)
    if (CMAKE_CXX_FLAGS MATCHES "/W[0-4]") # Force to always compile with W4
        STRING(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
elseif (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
endif()
