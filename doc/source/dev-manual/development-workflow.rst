====================
Development Workflow
====================

.. contents::

-----------
Description
-----------

The |Celeste| development cycle consists of "branch -> commit -> review -> merge" cycles for all future code changes.  To summarize:

1. The developer checks out a new branch
2. The developer commits code to the branch and pushes to the remote Source Control Management system (SCM)
3. The developer commits code to the branch
4. The developer opens a pull request (PR) for merging code changes in
5. Reviewers review the code and put comments/suggestions
6. The developer address fixes/issues/complaints/questions made by the code reviewers
7. Code reviewers approve the PR and merge the side branch in.

To maintain stable releases but allow for active development,two long-standing git branches are maintained:

* ``develop`` branch - This is where active development occurs
* ``master`` branch - This is where the release versions are maintained

Changes from the ``develop`` branch will be merged to ``master`` on a release-update basis.


----------------
Workflow Diagram
----------------

.. image:: workflow.svg


---------
Rationale
---------

Reasons for adopting this standardized workflow include, but is not limited:

1. Preventing developers from making bad changes to the main working branches
2. Unifying development efforts, as opposed to having developers branch off forever with different variants of the codebase for different fixes and features
3. Enforce knowledge transfer among the developers through code reviews, so that the project does not die
4. Bring some control and sanity to the development of |Celeste|
5. Enforce consistent quality in software delivery for |Celeste| users


-----------
Enforcement
-----------

To enforce this workflow, the following rules must be enabled in the hosting SCM:

* No direct pushes to ``master``, except by the integrator.
* No direct pushes to ``develop``, except by PR approval.
* PRs are blocked from the approval threshold unless the associated branch passes build and tests.
* The ``master`` and ``develop`` branches are locked from forced pushes / overwrites, but allows overwrites for all other branches.

Currently, |Celeste| is hosted on the public version of Bitbucket, which has limitations prventing us from enforcing the above rules.  There is an ongoing effort to (re)build the infrastructure using Jenkins and Gerrit in order to enforce these good software engineering practices.


-----------------------------
Workflow-Supporting Git Usage
-----------------------------

|Celeste| uses ``git`` for source code version control, and as such, developers should familiarize themselves with advanced ``git`` commands through :doc:`git` to support this workflow.


---------------------------------------------
Continuous Integration (CI) and Delivery (CD)
---------------------------------------------

While not in place at time of writing, Jenkins will be used for CI and CD.  Jenkins allows for the following to be done automatically in response to codebase updates:

1. Preventing PRs from being approved if the associated branch failed to build or failed tests
2. Enforcing coding style conformation
3. Automatically deploy new documentation and binaries for code updates to servers
