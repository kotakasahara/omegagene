# myPresto/omegagene

myPresto/omegagene is a Molecular Dynamics package (more description here).


## myPresto/omegagene Software Requirements

* CMake 3.4+
* For the GPU-version of myPresto/omegagene, CUDA 7.0+ is required for C++11 support.
* A C++11 compiler:
    * GCC: 4.8+
    * Clang: 3.6+
    * AppleClang: 5.0+
    * Intel: 15.0+ (minimal version required by CUDA 7.0+)

## Building myPresto/omegagene (Simple)

1. Set up a target build folder:

        # in <PROJECT_ROOT> directory
        $ mkdir target
        $ cd target

1. Configure the build.  CMake will determine all the external software dependencies for the selected build variant, and exit with errors if the dependency requirements are not met.  CMake must be invoked on the `CMakeLists.txt` file in the **<PROJECT_ROOT>** directory:

        # Run ONE of the following commands to configure for building the desired variant of myPresto/omegagene
        # in <PROJECT_ROOT>/target directory
        $ cmake -DCELESTE_WO_NS=1 ..
        $ cmake -DCELESTE_GPU=1 ..
        $ cmake -DCELESTE_GPUECP=1 ..
        $ cmake ..

1. Build the software:

        # The verbose flag is optional
	$ make VERBOSE=1

Additional build options and configurations can be found by invoking CMake `help`:

    # in <PROJECT_ROOT>/target
    $ cmake -LH ..
    ...
    ... help and custom build flag descriptions
    ...

## Documentation

Extensive documentation for celeste is available under the `/docs` directory.

To generate the documentation in HTML form:

    $ cd docs
    $ cmake .
    $ make
