Git help
********

This section gives the instructions of how to use the git from command line.
remember that you can use the git in the pycharm interface as well.

We assume that the 'origin' remote is the github repository. It is the default,
if it does not work read about git remote.

Issues
======

Issues are necessary to maintain a streamlined code when adding or modifying the codebase.

An issue is opened in the github interface, which gives it the number.
Then, open a branch with the name Issue<issue #> (issue # from git).

.. code-block::
    git branch -b Issue<#>

To return to the main version of the code:
.. code-block::
    git checkout master

And back to the issue:

.. code-block::
    git checkout Issue<#>


Updating the repository
=========================

Updating the code in  github is called pushing.

.. code-block::
    git push  origin


Updating the local code
=======================

Updating the local changes to the local repository
.. code-block::
    git commit -a

In order to get the changes from the github repository

.. code-block::
    git fetch Issue<#>
    git checkout Issue<#>

Merging the new master into your code
======================================

If the master changed during the development, it is wise to mege the current
branch with it. To do so, make sure you are in the Issue branch by typing

.. code-block::
    git branch

And make sure that the asterix is next to the Issue<#> branch.
Then,
.. code-block::
    git pull origin master

