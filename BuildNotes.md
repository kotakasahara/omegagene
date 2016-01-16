# Celeste Build Notes

Below is an assortment of build notes for handling different software dependencies and different platforms.

## General

### CUDA

* Library code that is built on top of CUDA must be built as **shared libraries**; otherwise, linker errors will appear during building on Linux platforms.  Hence, the entries in `CMakeLists.txt` specifying building CUDA-dependent libraries should be marked `SHARED` as such:

        CUDA_ADD_LIBRARY(CelesteFooCUDA SHARED foo.cu bar.cu)

* The version of gcc installed may be a later version than the the latest officially-supported host compiler for `nvcc`.  You will see an error like this:

        # in <PROJECT_ROOT>/target
        localhost:target local$ cmake -D CMAKE_C_COMPILER=/opt/local/bin/gcc-mp-5 -D CMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-5 -D CELESTE_GPU=1 ..
        ...
        ... generating Makefiles
        ...

        localhost:target local$ make
        ...
        ... building code
        ...
        /usr/local/cuda/include/host_config.h:115:2: error: #error -- unsupported GNU version! gcc versions later than 4.9 are not supported!

    While not recommended, this can be fixed by commenting out the appropriate `#error` macro in `<CUDA_ROOT>/include/host_config.h`:

        #if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 9)

        // #error -- unsupported GNU version! gcc 4.10 and up are not supported!

        #endif /* __GNUC__> 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 9) */


## Linux

### MPI

* Installation of MPICH or OpenMPI may not include adding `mpicc`/`mpic++` to the `$PATH`, resulting in the following error when `cmake` is invoked:

        Could NOT find MPI_C (missing: MPI_C_LIBRARIES MPI_C_INCLUDE_PATH)

    To resolve this, simply add the directory containing `mpicc`/`mpic++` to the `$PATH` in the `ENVIRONMENT` or in the `cmake` invocation:

        localhost:target local$ PATH=$PATH:/usr/lib64/mpich/bin cmake .....


## Mac OS X

* Unfortunately, with the newer versions of CUDA on the Mac, gcc is not a supported host compiler.  Attempting to compile CUDA code using gcc as the host compiler will result in an error message that looks like this:

        # in <PROJECT_ROOT>/target
        localhost:target local$ cmake -D CMAKE_C_COMPILER=/opt/local/bin/gcc-mp-5 -D CMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-5 -D CELESTE_GPU=1 ..
        ...
        ... generating Makefiles
        ...

        localhost:target local$ make
        ...
        ... building code
        ...
        nvcc fatal   : GNU C/C++ compiler is no longer supported as a host compiler on Mac OS X.

    Intel's ICC does not appear to be a supported host compiler for CUDA on Mac OS X either.  Only Clang appears to be a supported host compiler, but this only applies to "AppleClang" (the version of Clang maintained by Apple).  Attempting to use (newer versions of) mainline Clang will result in an error message that looks like this:

        # in <PROJECT_ROOT>/target
        localhost:target local$ cmake -D CMAKE_C_COMPILER=/opt/local/bin/clang-mp-3.7 -D CMAKE_CXX_COMPILER=/opt/local/bin/clang++-mp-3.7 -D CELESTE_GPU=1 ..
        ...
        ... generating Makefiles
        ...

        localhost:target local$ make
        ...
        ... building code
        ...
        nvcc fatal   : The version ('30700') of the host compiler ('clang') is not supported

    While not recommended for ABI/linking reasons, issues such as this above can be resolved by specifying a _different_ compiler as the host compiler for `nvcc`:

        # where /usr/bin/clang symlinks to AppleClang
        localhost:target local$ cmake -D CMAKE_C_COMPILER=/opt/local/bin/clang-mp-3.7 -D CMAKE_CXX_COMPILER=/opt/local/bin/clang++-mp-3.7 -D CELESTE_GPU=1 -D CUDA_HOST_COMPILER=/usr/bin/clang ..


## Windows

* MSVC does not define the alternative tokens for logical operators (i.e. `and` in place of `&&`) by default.  See http://stackoverflow.com/questions/24414124/why-does-vs-not-define-the-alternative-tokens-for-logical-operators.  This issue can be circumvented by including the following header in source files that use alternative tokens:

        #include <ciso646>

    The correct solution is to disable C++ language extensions in MSVC by use of the `/Za` compiler flag; however this flag is known to be buggy and will result in ODR errors during linking.  See the following articles:

    * http://cidebycide.blogspot.com/2015/10/visual-studio-2015-icu-and-error-lnk2005.html
    * http://stackoverflow.com/questions/31808256/multi-file-iostream-error-lnk2005-in-vs2015-with-za
