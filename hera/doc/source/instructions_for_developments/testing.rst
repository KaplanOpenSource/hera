Testing
*******

Testing Hera itself
===================

Smoke Test
----------

A Smoke test is a minimal test that things might be working. The ``smoke-test``
script in the ``hera`` main folder does exactly that:

- First, it checks to see that all the Python files in the project are,
  indeed, valid Python, by compiling them.

  Note that this does not even check that imports are correct, just
  syntactic validity.

- Then, it goes over the scripts in the ``./bin`` folder and try to execute
  each of them with the ``-h`` flag. This does invoke some import statements,
  but not much else

Internal Tests
--------------

Hera libraries include tests relying on the ``pytest`` library and toolset.
The way to run these tests is::

   $ pytest

The tool detects the tests automatically and runs them.

A more comprehensive test also executes all the Jupyter Notebooks in the
documentation folders. This is only available if you install the documentation
requirements -- i.e. you need to have run::

   $ pip install -r requirements-doc.txt

With these in place, you can invoke pytest with the flags::

   $ pytest --nbmake

At the time of writing, and on the computer used to write this,
this takes about 3:30 minutes. However, it produces a lot of output,
which may be confusing. A more minimalist variation is the command::

   $ pytest --nbmake --nbmake-find-import-errors

The ``--nbmake-find-import-errors`` flag limits the errors being reported,
leaving only import errors, which typically also have shorter tracebacks.