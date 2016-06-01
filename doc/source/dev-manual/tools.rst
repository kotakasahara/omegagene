====================================
Tools Used for |Celeste| Development
====================================

.. toctree::
   :maxdepth: 2

   doxygen
   jenkins

----------------------
Source Code Management
----------------------

**git**
  |Celeste| uses git as the version control system.  Currently, it is
  hosted on `Bitbucket <https://bitbucket.org/kotakasahara/celeste>`_.
  Instructions for setting up git for |Celeste| can be found in the project
  README Other basic tutorial material for ``git`` can be found on the web.

**Gerrit [FUTURE]**
  All code changes should go through a code review system.  Currently this
  role is fulfilled by BitBucket, but BitBucket is not flexible enough to handle
  stringent code reviews.

**Jenkins [FUTURE]**
  All changes pushed to Gerrit are automatically compiled and otherwise
  checked on various platforms using a continuous integration system
  [to be determined].

**Redmine [FUTURE]**
  Bugs and issues, as well as some random features and discussions,
  are tracked at http://redmine.gromacs.org.


------------
Build System
------------

**CMake**
  CMake is our main tool used in the build system.

**CPack**
  [FUTURE] Main tool used to package installers

**CTest**
  [FUTURE] Main tool used to package installers

**Docker**
  [FUTURE] Use a docker host to dynamically provision a slave, run a single
  build+test in isolation, then tear-down that slave.


----------------------------
Code Analysis and Formatting
----------------------------

**cppcheck**
  [FUTURE]

**clang static analyzer**
  [FUTURE]

**clang-format**
  [WRITEME]

**include-what-you-use**
  [FUTURE]


------------------------
Documentation Generation
------------------------

Building the |Celeste| documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. TODO: Move this onto a separate page

For now, there are multiple components, formats and tools for the
|Celeste| documentation, which is aimed primarily at version-specific
deployment of the complete documentation on the website.

This is quite complex, because the dependencies for building the
documentation must not get in the way of building the code
(particularly when cross-compiling), and yet the code must build and
run in order for some documentation to be generated. Also, man page
documentation (and command-line completions) must be built from the
wrapper binary, in order to be bundled into the tarball.

The outputs of interest to most developers are generally produced in the
``docs/html/`` subdirectory of the build tree.

The following make targets are the most useful:

``manpage``
  Makes man pages from the RST files with Sphinx
``html``
  Makes html pages from the RST files with Sphinx
``doxygen``
  Makes the code documentation with Doxygen


**Doxygen**
  `Doxygen <http://www.doxygen.org>`_ is used to extract documentation from
  source code comments.

**graphviz (dot)**
  The Doxygen documentation uses ``dot`` from `graphviz
  <http://www.graphviz.org>`_ for building some graphs.  The tool is not
  mandatory, but the Doxygen build will produce warnings if it is not
  available, and the graphs are omitted from the documentation.

**mscgen**
  The Doxygen documentation uses `mscgen
  <http://www.mcternan.me.uk/mscgen/>`_ for building some graphs.  As with
  ``dot``, the tool is not mandatory, but not having it available will result
  in warnings and missing graphs.

**Doxygen issue checker**
  [FUTURE]

**module dependency graphs**
  [FUTURE]

**Sphinx**
  `Sphinx <http://sphinx-doc.org/>`_; at least version 1.2.3) is used
  for building some parts of the documentation from reStructuredText
  source files.

**LaTeX**
  Also requires ImageMagick for converting graphics file formats.

**linkchecker**
  [FUTURE]