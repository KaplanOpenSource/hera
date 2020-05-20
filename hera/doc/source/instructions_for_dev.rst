Instructions for developers
###########################


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

In python 3.6 it is possible to format a string with the 'f""' directive:

.. code-block:: python

    a = 5
    str = f"The value of a is {a}"


will result in str= The value of a is 5


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


