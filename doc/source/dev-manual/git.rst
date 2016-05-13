=========================
Select Advanced Git Usage
=========================

While widely used, ``git`` is not always well understood and often is used by
developers as one would use ``svn``.  Below are some advanced usages of ``git``
that |Celeste| developers should be familiarize themselves with.  More
documentation is always available at the
`official Git Documentation <https://git-scm.com/docs>`_.

-------------------------------
Simple Rebasing: ``git rebase``
-------------------------------

Git allows developers to rebase branches.  To visualize:

Before:

.. code-block:: text

        A---B---C topic
       /
  D---E---F---G master


Running ``git rebase``:

.. code-block:: bash

  git rebase master topic


After:

.. code-block:: text

                  A'--B'--C' topic
                 /
    D---E---F---G master


**Developers should be rebasing their changes prior to opening a PR.**  This
will help prevent merge conflicts and the resulting non-linear history that
appears when PR is merged in, which makes future code inspection of commits
unnecessarily difficult:

.. code-block:: text

  # BAD
        A---B---C
       /         \
  D---E---F---G---H master

  # GOOD
  D---E---F---G---A'--B'--C' master


If the commits in the side branch were already pushed prior to rebase, they
will need to be force-pushed after the rebase.

---------------------------------------
Interactive rebasing: ``git rebase -i``
---------------------------------------

In addition, ``git`` supports interactive rebasing.  This allows develoeprs to:

1. Squash related commits together into one commit
2. Edit commit messages
3. Delete non-useful commits

Before:

.. code-block:: text

        A---B---C---D----E---F topic
       /
  D---E---F---G master


Running ``git rebase -i master topic`` will open an interactive ``vim`` shell.
Simply edit the file as instructed and save to disk to execute the rebase
 operations:

.. code-block:: bash

  pick 358bc5a commit A in topic branch
  squash 1b01e0a commit B in topic branch
  pick 755189f commit C in topic branch
  pick b95b8a4 commit D in topic branch
  squash 084bf71 commit E in topic branch
  squash 01e4fde commit F in topic branch

  # Rebase 897f795..084bf71 onto 897f795 (5 command(s))
  #
  # Commands:
  # p, pick = use commit
  # r, reword = use commit, but edit the commit message
  # e, edit = use commit, but stop for amending
  # s, squash = use commit, but meld into previous commit
  # f, fixup = like "squash", but discard this commit's log message
  # x, exec = run command (the rest of the line) using shell
  #
  # These lines can be re-ordered; they are executed from top to bottom.
  #
  # If you remove a line here THAT COMMIT WILL BE LOST.
  #
  # However, if you remove everything, the rebase will be aborted.
  #
  # Note that empty commits are commented out


After:

.. code-block:: text

        (AB)---C---(DEF) topic
       /
  D---E---F---G master


-------------------------------------------------------------------
Applying Exising Commits to a Different Branch: ``git cherry-pick``
-------------------------------------------------------------------

Sometimes, we want to apply a desirable change from one working branch to
another working branch.  In the below example, would like to apply commit ``C``
to the ``master`` branch.

Before:

.. code-block:: text

        76cada---62ecb3---b886a0 topic
       /
  dd2e86---946992---9143a9---a6fd86---5a6057 master


Running ``git cherry-pick``:

.. code-block:: bash

  git checkout master
  git cherry-pick 62ecb3


After:

.. code-block:: text

        76cada---62ecb3---b886a0 topic
       /
  dd2e86---946992---9143a9---a6fd86---5a6057--[62ecb3] master
