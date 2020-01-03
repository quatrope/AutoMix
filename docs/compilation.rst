.. _compile:

Compilation
-----------

To compile the AutoMix sampler for the examples included within the AutoMix package, the user should edit the `Makefile` (supplied) so that the appropriate C compiler is used (currently set to be the GNU sampler gcc). Other aspects of the Makefile may be edited if required (for example the compile flags, or libraries) 

Optimized Compilation
^^^^^^^^^^^^^^^^^^^^^

Typing::

    make all

in the shell, in the AutoMix folder where the AutoMix distribution was unzipped to, will compile all the programs that were distributed with this package.

Compile with Debug Information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To compile all programs with debug information (`-g` flag enabled), type::

    make alldebug

instead. The supplied programs have all been debugged but we acknowledge that use of a debugger can often help to understand how the program works. 

Any of the programs can also be made individually, with the appropriate make command.

Example
^^^^^^^

For example::

    make amtoy1

will make amtoy1. See the Makefile for further examples.

The executables have the form

* `amNAME` (for the optimised version)
* `amNAMEd` (for the version that can be debugged)

where `NAME` is the name of the example.

Clean
^^^^^

To remove the executables and object (.o) files, type::

    make clean   

To compile the sampler for a new example, the Makefile can be edited to contain the appropriate commands (remembering to include command for compiling any dependencies) for the new program.
