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

Hera libraries include tests relying on the `pytest`_ library and toolset.
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

.. _`pytest`: https://pytest.org

Writing Tests
-------------

As mentioned above, the tests in Hera are based on `pytest`_. It is highly
recommended to read some of the explanations and how-ro guides in its
documentation before writing tests, but we will include a brief introduction
here.

`pytest`_ uses a few simple conventions, which we'll usually follow, to make things
easier.

What is a test
..............
Any module (file) whose name begins with ``test_`` is considered a test
module, and in a test module, any function whose name begins with ``test_``
is a test.

What is inside a test
.....................
As the pytest documents `explain`_, a test is usually done by invoking some
action and then making assertions about its result; but this sometimes also
requires preparations and cleanups.

The actions are, basically, your code.

For assertions, we mostly rely on the Python ``assert`` statement, which
checks a condition. Thus, if we want a test which verifies that basic
arithmetic works, we might write::

    def test_two_plus_two():
        result = 2 + 2      # Action
        assert result == 4  # Assertion

Preparations and cleanups have two qualities which we can use to improve
our testing code.

First, it is typical for the same set of preparations and cleanups to
be used in many tests; in these cases, we'd want to factor them out of
the tests so they don't have to be repeated.

Second, usually preparations and cleanups are related -- that is, the
preparation sets up some objects for the test to use (e.g. some file
to be read), and the cleanup needs to tear down the same objects (remove
the file). It would be nice if we could somehow bundle them together,
so that whenever the preparation is done, the cleanup will be done as
well.

In pytest, the common way to handle preparations and cleanup is using
`fixtures`_. They are a complex and powerful construct, and the document
in the link is long; but the essentials are:

- You define a fixture by defining a function and decorating it with
  ``@pytest.fixture``
- You use a fixture in a test by giving its name as an argument to
  the test function. When the test runs, the return value from the
  fixture function will be supplied as an argument
- If your preparation should imply a cleanup, use ``yield`` instead of
  ``return``, and write the cleanup part following the ``yield`` (as
  shown `here`_).

Fixtures can be used for more than that, and we do -- see in the next
section.

.. _explain: https://docs.pytest.org/en/latest/explanation/anatomy.html
.. _fixtures: https://docs.pytest.org/en/latest/how-to/fixtures.html
.. _here: https://docs.pytest.org/en/latest/how-to/fixtures.html#yield-fixtures-recommended

Special pytest features for Hera
................................

For work with Hera, we have set up several special tools and behaviors.

Test Database
,,,,,,,,,,,,,

When you run tests using ``pytest``, the tests can use mongodb through
Hera's toolkits and APIs. However, they will not be running against your
regular database; they run against a special database that is set up
just for the tests.

This allows you to use real database-modifying operations in your tests,
without fear of messing with your production data.

Fixtures for regression tests
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
We include a pytest plug-in for regression tests -- tests whose action
generates some data, and the assertion is that the data hasn't changed
since the last time we ran the test. This is based on `pytest-regressions`_
with a small addition to make it more comfortable to use with data as
stored in Mongo; our version is called ``cached_data_regression``.

As the subtitle says, the idea is a fixture -- an object which the test
receives as an argument, and has a method to check the data. Behind
the scenes, the tool saves the data in a file, and will read it from
there (to compare) the next time.

For this to work, the first time you run such a test it will fail (because
the data isn't there yet), but will write the data; and you can ask pytest
explicitly to regenerate by giving the flag ``--regen-all``.

.. _pytest-regressions: https://pytest-regressions.readthedocs.io/en/latest/overview.html

See it all come together
........................
See an example for tests written with regression-test fixtures, and also
with module-specific fixtures, in ``hera/measurements/GIS/test_shapes.py``.

In that file:

- Two fixtures are defined, one to fetch the shapes toolkit, and the other
  to select a shape file to work with

- And then three tests, each using the two fixtures:

  + One uses a simple assertion
  + One uses geopandas' ``to_json()`` method to format data in a way the
    "vanilla" pytest-regressions can handle
  + And the last one uses our special ``cached_data_regression`` fixture

The data files are stored in a folder next to the test module,
``hera/measurements/GIS/test_shapes``, with one file for each test.
