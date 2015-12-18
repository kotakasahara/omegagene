# Celeste

Celeste is a Molecular Dynamics package (more description here).


## Celeste Software Requirements

* A C++11 compiler (clang 3.2+, gcc 4.9+)
* CMake 3.2+
* For the GPU-version of Celeste, CUDA 5.0+ is required


## Building Celeste
1. Set up a target build folder

        # from project root directory
        mkdir target
        cd target

1. Configure the build

        # Run ONE of the following commands to configure for building the variant of celeste that you want
        # in /target
        cmake -DCELESTE_WO_NS=1 ..
        cmake -DCELESTE_GPU=1 ..
        cmake -DCELESTE_GPUECP=1 ..
        cmake ..

1. Build the software

        # The verbose flag is optional
        make VERBOSE=1

Additional build options and configurations can be found by invoking CMake help:

    # in /target
    cmake -LH ..
    ....
    ....help and custom flag descriptions....
    ....
