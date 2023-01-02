Technical Support
#################


Synchronize intermediate version
================================

In order to synchronize the code with another user
just type::

git pull origin [developer name]

In order to avoid conflicts, you can either accept their version::

git pull origin [developer name] -X theirs

or force your own version::

git pull origin [developer name] -X ours



Python
=======

Formatting string
-----------------

In python >= 3.6 it is possible to format a string with the 'f""' directive:

.. code-block:: python

    a = 5
    s = f"The value of a is {a}"


will result in ``s`` holding ``"The value of a is 5"``

Python 3.8 added the ``=`` suffix for debug printings:

.. code-block:: python

    a = 5
    s = f"Now {a=}"

will result in an error with Python < 3.8, but newer
Pythons will give ``"Now a=5"``

Debug logging
-------------

To enable debug in specific loggers, see
:py:func:`hera.utils.logging.helpers.initialize_logging`.

.. automodule:: hera.utils.logging.helpers


Creating objects in real-time
-----------------------------

Sometimes (very rarely) it is important to derive objects in realtime.
In Hera we use it to create object for the ORM (object relation mapping) of MongoDB.

For example to create a new class that derives form myFather:

.. code-block:: python

    newClass = type('mynewclass', (Son,Father,), {})

    newClassInstance = newClass()

is Equivalent to

.. code-block:: python

    class mynewclass(Son,Father):
        ...


