================================
|Celeste| Infrastructure Roadmap
================================

.. contents::

The |Celeste| codebase and infrastructure needs to grow and advance far beyond
its current capabilities to be able to survive as a project and be useful
outside of a group of projects under one or two research groups.  The roadmap
focuses on four broad areas of development.

--------------------------
Project Tooling / Pipeline
--------------------------

^^^^^^^^^^^^^^^^^^^^^^^^^^
Expanding the CMake System
^^^^^^^^^^^^^^^^^^^^^^^^^^

While the project's current build system is CMake, more work needs to be done
to ensure production-level quality and higher levels of build automation.
This includes:

* Softcode all software versioning / build information into CMake variables
  across both source and docs (Sphinx and Doxygen)
* More robust CMake scripts to support more OS platforms and future build
  configurations
* CPack integration, so that we can build release tarballs/installers

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Automated Codebase Checks/Cleaning
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following tools need to be integrated into the build system as custom build
targets for automated code checks/cleanups (in no particular order of priority):

* ``clang-format``/``uncrustify`` - for enforcing code formatting styles
* ``clang-tidy`` - for checking common C++ style/pattern errors
* ``clang-static-analyzer``/``cppcheck`` - for running static analysis checks
  on the codebase

^^^^^^^
Doxygen
^^^^^^^

Currently there are neither docstrings in the codebase nor Doxygen integration
to take advantage of them.  Doxygen integration needs to be added to the build
system.

^^^^^^^^^^^^^^^^
Jenkins + Docker
^^^^^^^^^^^^^^^^

Developers should not be concerned with manually testing to approve PRs, nor
should they be concerned with the software delivery process.  Jenkins + Docker


^^^^^^
Gerrit
^^^^^^

We need a code review system that actually enforces the workflow policies
documented in :doc:`documentation-workflow`.


-----------------
Language Upgrades
-----------------

^^^^^^^^
C++11/14
^^^^^^^^

The current codebase is in C++98.  We should try to move the codebase to C++11, if not C++14.  C++11/14 will provide the following:

1. simplify the codebase (shorter code)
2. a lot more code optimizations possible than C++98.

This needs to be done carefully, since some compilers like Intel icpc do not fully support C++11/14 yet.  I am also not sure what versions of gcc/clang are on the HPC clusters celeste is running on, nor the difficulty in upgrading compiler versions - we need gcc >= 5.1 or clang >= 3.4 (3.7 if you want openmp).


^^^^^^^^^^^^^^^^^^^^^
Remove C++ non-idioms
^^^^^^^^^^^^^^^^^^^^^

The current codebase is contains many Fortran and C-style programming that prevent the codebase from taking advantage of C++ efficiency.  Examples include but are not limited to:

1. Manual memory allocation/deallocation: RAII is ignored, which prevents compiler optimizations.  Using RAII will significantly reduce the current codebase size, since a large portion of the code goes into writing constructors and destructors.

2. Java Object model in C++: For example, from what I can tell, the main purpose of CelesteObject is to store constants.  Whenever possible, we should try to avoid parent objects and inheritance/virtual-functions because they are slow, and we should be using namespaces and templates instead.

3. The "triple" type in CelesteTypes.h: While I don't see it actually being used in the codebase, I think this is deprecated by the C++11 "tuple" type, which allows for *any* number of types in the struct.

The hope is that celeste will be much faster when the code is changed to proper C++-style programming.


-----------------------
Codebase Modularization
-----------------------


This is a long-term "high-level" goal/idea (Please let me know if this is not a good idea). I think the issue is two-fold:

1. the current codebase is monolithic (no modules/namespacing etc), and so new features are hard to add/debug
2. the style of MD development has been for each researcher to (re)build components.  For example, in tplgene, there are at least 3 source files doing PDB-parsing.

From an MD developer's point of view, it would be nice if celeste provided a collection of libraries that one can choose from to easily build his own custom MD "frontend"/algorithm/workflow.  This is similar to the LLVM + clang compiler strategy.  I believe OpenMM (Vijay Pande's group) also takes this approach, and I think the celeste python toolkit is striving for that direction (at the scripting level) as well : ).

I've discussed this with Gert-Jan before as to what needs to be done:

1. Build standard libraries for I/O and standardized wrappers for talking to 3rd party libraries
2. Standardize the data structures so that different algorithms can be easily interchanged (PDF documentation).
3. Build standard libraries for the most common MD analysis codes (i.e. PCA)
4. (If possible) C bindings for the libraries so that users can directly use the simulation code from python

I think the benefits for partitioning celeste into component libraries are:
1. celeste code development can be done in parallel
2. it will be much easier to maintain/document the code
3. it will be much easier for newcomers to learn and add improvements to the code.


As an example, the sister project at https://bitbucket.org/bjma/noname has some "standard" libraries for efficiently reading and writing PDB and TPL files.

Since celeste is still early in the stages of development, I think building libraries is a good move in that it will allow celeste to evolve faster and encourage external more developers to join and contribute.  It will become much more difficult to do so when celeste grows to an even larger codebase than it currently is.


------------------
External Libraries
------------------

^^^^
HDF5
^^^^

HDF5 is a standard file format for large scientific/numerical data.  I think using this format as the MD trajectory/output as opposed to other/custom formats might be helpful for the following reasons:

1. HDF5 is an HPC industry standard, and can support multi-terabyte-sized files
2. HDF5 data is structured
3. HDF5 can store different float/int sizes, arrays, strings, jpeg images, etc - it is a filesystem-in-a-file, which is useful if you want to store every information about an MD run in one file
4. HDF5 has its own efficient data compression mechanisms
5. There are libraries for opening large HDF5 files in parallel through MPI - this might be useful for post-MD analysis codes.

^^^^
HDF5
^^^^
