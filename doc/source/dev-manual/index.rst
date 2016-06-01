=========================
Celeste Developers' Guide
=========================

The developers' documentation consists of two parts:

* A developer guide that provides an overview of the |Celeste| codebase, and
  includes more detailed resouces such as coding guidelines, infrastruture
  setup, and information on tools used during development.
* Doxygen documentation extracted from comments in C/C++ code, documenting the
  actual C/C++ code.

The developer guide focuses on things that are tightly tied to the code
itself, such as helper scripts that reside in the source repository and
organization of the code itself, and requires the documentation to be
updated in sync.  The guide is currently split into a few focus areas:

* Overview of the |Celeste| codebase
* Overview of some important implementation aspects.
* Overview of supporting infrastructures, such as continuous integration
* Guideline of the development process, and what workflows are to be expected
* Guidelines to follow when developing |Celeste|
* Summary of tools used, and how to use them.

The documentation does not yet cover all areas of |Celeste| development, but
more content is being (slowly) added.

Contents:

.. toctree::
    :numbered:
    :maxdepth: 2

    tools
    development-workflow
    style
    git
    project-roadmap
    infrastructure-roadmap
..    build-system
