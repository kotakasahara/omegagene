==========================
Code Formatting Guidelines
==========================

Enforcing a consistent formatting has a few advantages:

* No one needs to manually review code for most of these formatting issues,
  and people can focus on content.
* A separate automatic script (see below) can be applied to re-establish the
  formatting after refactoring like renaming symbols or changing some
  parameters, without needing to manually do it all.

Below are the guidelines for formatting code in the various languages used
for |Celeste| development.

-------------------------
C++ Formatting Guidelines
-------------------------

The following list provides the general formatting/indentation rules for
|Celeste| code (C/C++):

* Basic indentation is four spaces.
* Keep lines at a reasonable length.  Due to the natural complexity of the code,
  use 120 characters as a guideline.  If you end up indenting very deeply,
  consider splitting the code into functions.
* Never use tabs; use only spaces.  Most editors can be configured to generate
  spaces even when pressing tab.  Tabs (in particular when mixed with spaces)
  easily break indentation in contexts where settings are not exactly equal
  (e.g., in ``git diff`` output).
* No trailing whitespace.
* ``extern "C"`` and ``namespace`` blocks are not indented, but all others
  (including ``class`` and ``switch`` bodies) are.

Most of the above guidelines are enforced using ``clang-format``, an automatic source
code formatting tool.  While it is a suitable tool for enforcing the basic rules,
``clang-format`` may not have enough formatting options to serve our needs over time,
and so there is consideration to migrating over to ``uncrustify``, which provides
over 600 formatting rules.


----------------------------
Python Formatting Guidelines
----------------------------

Please refer to the formatting rules in the Python official
`PEP 8 Guidelines <https://www.python.org/dev/peps/pep-0008>`_.


-------------------------------
Copyright Formatting Guidelines
-------------------------------

[WRITEME]

.. Additionally:

.. * All source files and other non-trivial scripts should contain a copyright
..   header with a predetermined format and license information (check existing
..  files).  Copyright holder should be "the Celeste development team" for the
..  years where the code has been in the |Celeste| source repository, but earlier
..  years can hold other copyrights.
.. * Whenever you update a file, you should check that the current year is listed
..  as a copyright year.

.. The copyright guidelines are enforced by a separate
.. Python script.  See :doc:`uncrustify` for details.  Note that due to the
.. nature of uncrustify (it only does all-or-nothing formatting), it enforces
.. several additional formatting rules in addition to those above.
