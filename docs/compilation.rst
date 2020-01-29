.. _compile:

Compilation
===========

Library Compilation
-------------------

To compile the automix library, simply type ``make`` on the command line::

    $ make

This default command compiles only the library with optimization flags ``-O3`` on.

If you need to link your program against a library with debug information, type::

    $ make DEBUG=1

instead.

If the compilation went without problems, you should see the ``libautomix.so`` file created.

Typically, libraries are installed in the ``/usr/local/lib`` directory,
and include files such as ``automix.h`` are installed in ``/usr/local/include``.
If the default install location is enough for your needs, type::

    $ sudo make install

``sudo`` privilege may be needed to copy files into the ``/usr/local`` directory.

If you prefer to install in a different location you can provide it with the ``PREFIX`` command line flag::

    $ make install PREFIX=/my/preferred/path/

Compiling your program against libautomix
-----------------------------------------

If your program is contained in a single ``main.c`` file, and assuming no extra libraries or include files are required::

    $ cc main.c -lautomix -o myprogram

should be enough to compile. If ``/usr/local/`` is not a default search place for your libraries,
add a ``-L/usr/local/lib`` and a ``-I/usr/local/include`` to the compilation line::

    $ cc main.c -L/usr/local/lib -I/usr/local/include -lautomix -o myprogram

You should have a ``myprogram`` executable compiled.

Compiling test and tutorial
---------------------------

Along with the library you can compile and run the tests::

    $ make test
    $ ./test

Or the tutorial provided in this documentation::

    $ make tutorial
    $ ./tutorial

Compiling "user*" legacy programs
---------------------------------

.. warning::

  The user examples are deprecated on version 2.x.
  If you want to compile AutoMix using the old interface with ``user*`` files, please refer to version 1.x of this package.

The original AutoMix was program-oriented and offered several example programs.
To compile the example programs type::

    $ make examples

After compilation you should see the programs
``amtoy1``, ``amtoy2``, ``amcpt``, ``amcptrs``, ``amddi``, and ``amrb9``.

You can compile each one of them individually as well, for example::

    $ make amtoy1
    $ ./amtoy1

Clean
-----

To remove the executables and object (.o) files, type::

    make clean   
