# Celeste

Celeste is a Molecular Dynamics package (more description here).


## Celeste Software Requirements

* CMake 3.2+
* For the GPU-version of Celeste, CUDA 7.0+ is required for C++11 support.
* A C++11 compiler:
    * GCC: 4.8+
    * Clang: 3.6+
    * AppleClang: 5.0+
    * Intel: 15.0+ (minimal version required by CUDA 7.0+)
* OpenMP 3.1+


## Building Celeste

This is non-exhaustive description for building Celeste.  For platform-specific details, please refer to the [Build Notes](BuildNotes.md).

1. Set up a target build folder:

        # in <PROJECT_ROOT> directory
        localhost:celeste local$ mkdir target
        localhost:celeste local$ cd target

1. Configure the build.  CMake will determine all the external software dependencies for the selected build variant, and exit with errors if the dependency requirements are not met.  CMake must be invoked on the `CMakeLists.txt` file in the **<PROJECT_ROOT>** directory:

        # Run ONE of the following commands to configure for building the desired variant of celeste
        # in <PROJECT_ROOT>/target directory
        localhost:target local$ cmake -DCELESTE_WO_NS=1 ..
        localhost:target local$ cmake -DCELESTE_GPU=1 ..
        localhost:target local$ cmake -DCELESTE_GPUECP=1 ..
        localhost:target local$ cmake ..

1. Build the software:

        # The verbose flag is optional
        localhost:target local$ make VERBOSE=1

Additional build options and configurations can be found by invoking CMake `help`:

    # in <PROJECT_ROOT>/target
    localhost:target local$ cmake -LH ..
    ...
    ... help and custom build flag descriptions
    ...
