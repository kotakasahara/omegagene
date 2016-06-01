================================
|Celeste| Infrastructure Roadmap
================================

.. contents::

The |Celeste| codebase and infrastructure needs to grow and advance far beyond
its current capabilities to be able to survive as a project and be useful
outside of a group of projects under one or two research groups.  The roadmap
focuses on five broad areas of development.


-----------------------
Project Support Tooling
-----------------------

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
* CPack integration to easily build release tarballs/installers

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
to take advantage of them.  Doxygen needs to be integrated into the build
system to automatically generate documentation from docstrings to be added.


------------------------------------------
Continuous Integration / Delivery Pipeline
------------------------------------------

^^^^^^^^^^
Unit Tests
^^^^^^^^^^

As the |Celeste| codebase grows, there is no scalable method for verifying that
code changes don't break existing functionality other than writing tests.  We
should be using `GTest <https://github.com/google/googletest>`_ and CTest to
facilitate automated testing of the codebase.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Automated Software Packaging
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Work needs to be done to automate the tarball/installer building process using
CPack, which comes with CMake.

^^^^^^^^^^^^^^^^
Jenkins + Docker
^^^^^^^^^^^^^^^^

Jenkins is a continuous integration (CI) tool that listens to git events
(i.e. pushes, pull requests openings, merges) and execute pre-determined scripts
based on those events.  Using all the aforementioned automation tools,
Jenkins can:

* Automatically build and test the codebase for regressions
* Automatically run code formatting / sanitizer / static analysis / other
  codebase checks (in parallel)
* Automatically update documentation server pages when a release is made
* Automatically update software tarballs/installers on the server when a release
  is made

Developers should never be concerned with manually testing to approve PRs, nor
should they be concerned with the software delivery process, and Jenkins will
address exactly those needs.

^^^^^^
Gerrit
^^^^^^

We need a code review system that works in conjunction with Jenkins and actually
enforces the workflow policies documented in :doc:`development-workflow`.


-----------------
Language Upgrades
-----------------

^^^^^^^^
C++11/14
^^^^^^^^

The current codebase is C++98-compatible but is compiled with the ``--std=c++11``
flag.  The intent is to break compatibility and be at least C++11-compatible,
if not C++14-compatible.  In addition to providing more expressive language
features, C++11/14 will bring the following benefits to |Celeste|:

* Simplification of the codebase (less/shorter code)
* Enabling of a lot more code optimizations possible than what is allowed
  in C++98.

Modernizing the codebase needs to be done carefully, since we will be running
into the following issues:

* Some compilers like Intel ``icpc`` do not fully support C++11/14 yet.
* Some compilers provide C++11/14 support in a version that is incompatible
  with existing libraries/tools.  For example, ``gcc`` 5.1 is not supported by
  CUDA 7.5.
* HPC clusters may not have the required compiler versions that are needed to
  fully compile C++11/14.  This cannot be avoided by simply installing the
  needed compiler to the user's home directory, because the code will be running
  on compute nodes, not the HPC cluster's login node.


^^^^^^^^^^^^^^^^^^^^^^
Remove C++ anti-idioms
^^^^^^^^^^^^^^^^^^^^^^

The current codebase is contains many Fortran and C-style programming idioms
that prevent the codebase from taking advantage of C++ efficiency.  Examples
include but are not limited to:

1. **Manual Memory Allocation/Deallocation**: C++ RAII is ignored, which
prevents compiler optimizations.  Using RAII will significantly reduce the
current codebase size, since a large portion of the code goes into writing
constructors and destructors.

2. **Java-style OOP Design**: Namespaces and templates are preferred over class
hierarchies for storing constants (i.e. ``CelesteObject``) and method
signatures.  Avoid parent objects and inheritance/virtual-functions because
they are slow.

3. **Not using STL constructs**: For example, ``CelesteTypes.h`` defined the
``triple`` type when this construct is available in C++11 in the more generic
``std::tuple`` type.

The hope is that |Celeste| will be much faster when the code is changed to
proper C++-style programming.


-----------------------
Codebase Modularization
-----------------------

The current issue is two-fold:

1. The current codebase is monolithic (no modules/namespacing etc), and so new
   features may be hard to add and/or debug.
2. Up until now, the style of MD development has been for each researcher to
   duplicate and/or (re)build components.  For example, in ``tplgene``, there
   are at least 3 source files/modules that perform the same PDB-parsing
   routines.

For the above reasons, development will not be scalable and organic, and it
would be beneficial if |Celeste| were developed instead as a collection of
libraries that developers can mix and match from to easily build their own
custom MD "frontend"/algorithm/workflow.  This is the architectural strategy
that the LLVM + clang compiler project follows.
The `OpenMM Project <http://openmm.org>`_ by Vijay Pande's research group also
takes this approach. Modularizing |Celeste| into componenet libraries will
provide for the following benefits:

* |Celeste| code development can be done in parallel.
* Maintanence and documentation efforts will be scalable.
* It will be much easier for newcomer developers  to learn and add
  improvements to the code.

Benson has discussed this with Gert-Jan before as to what needs to be done:

* Build standard libraries for I/O (i.e. reading and writing PDB/TPL files) and
  standardized wrappers for talking to 3rd-party libraries
* Standardize data structures so that different algorithms can be easily
  interchanged
* Build standard libraries for the most common MD analysis codes (i.e. PCA)
* C bindings for some of the libraries so that users can directly use the
  code from Python


--------------
New Directions
--------------

^^^^^^^^^^^^^
HDF5 Adoption
^^^^^^^^^^^^^

HDF5 is a standard file format for storing large amounts of scientific/numerical
data.  This format should be adopted by |Celeste| as the MD trajectory/output
as opposed to other/custom formats for the following reasons:

1. HDF5 is an HPC industry standard, and can support multi-terabyte-sized files.
2. HDF5 data is structured.
3. HDF5 can store different float/int sizes, arrays, strings, jpeg images, etc -
   it is a filesystem-in-a-file, which is useful for storing every available
   information about an MD run in one file.
4. The HDF5 library come with its own efficient data compression mechanisms.
5. There are libraries enabling parallel I/O for HDF5 files.  This will be
   useful for building post-MD analysis pipelining.
