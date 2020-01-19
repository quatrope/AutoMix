.. _run:

Running AutoMix
---------------

The sampler is run by typing the name of the executable, followed by run-time flags separated by spaces.
The run-time flags control a number of options within the sampler.
If flags are not supplied default values are used.
The flags can be used in any order.

Command-line Arguments
^^^^^^^^^^^^^^^^^^^^^^

The flags can be summarised as follows (I is assumed to be a positive integer):

.. option:: -m <mode>

   Controls the mode of the sampler.
     - Mode 0 is mixture fitting.
     - Mode 1 skips stage 1 and 2 if a file containing the mixture parameters is supplied.
     - Mode 2 fits AutoMix version of AutoRJ sampler (see Green, 2003 - full reference in thesis).

     Default uses `D=0`.

.. option:: -n <iters>

   Run the sampler for max(iters, nkk * 10000, 100000) iterations in the stage 1 RWM for each model k.

   Default value is 100,000.

.. option:: -N <iters>

    Run the sampler for `iters` Reversible jump iterations in stage 3.

    Default uses I=100000.

.. option:: -s <seed>

    Initialises the random number generator with seed I.

    Default uses clock as seed.

.. option:: -a <adapt>

    Controls whether or not adaptation is done in stage 3 RJ.

    If `adapt=0` no adaptation is done, if `adapt=1` adaptation is done.
    Default is 1.

.. option:: -p <perm>

    Controls whether or not random permutation is done in stage 3 RJ.

    If `perm=0` no permutation is done, if `perm=1` permutation is done.
    Default is 0.

.. option:: -t <dof>

    Controls whether standard Normal or T distributed variables are used in RWM and in RJ moves.

    If `dof=0` Normal variables are used, otherwise t-distributed variables with dof degrees of freedom are used.
    Default is 0.

.. option:: -f <basestring>

    Uses the string `basestring` as the basis for filenames (e.g. if `basestring=output`, filenames will be named `output_log.data`, `output_mix.data` etc).

    Default is "output".

.. option:: [-h, --help]

    Prints help information on command line arguments and exits.

.. option:: [-v, --version]

    Prints AutoMix version number and exits.

Example
^^^^^^^

As an example, typing::

    amtoy1 -m 0 -N 1000000 -p 1 -f toy1

runs the optimized mixture fitting version of the toy1 problem (see thesis, section 5.5.1) with 1 million RJ sweeps, enabling permutation and storing the output in files of the type `toy1_*.data`.
Running the sampler produces a summary of how the run is progressing.

For each of the models:

* In **stage 1** a countdown of the number of iterations remaining is printed to screen;
* In **stage 2** a summary of the mixture fitting is printed to screen.
  This summary consists of a countdown of the number of components in the current mixture, with the iteration number that the last component 
  was removed and an indicator `n` if the component was annihilated naturally, and `f` if the annihilation was forced.
* In the RJ **stage 3** a countdown of the number of iterations remaining is printed to screen. 
  No summary statistics are printed to screen. Instead all output from the sampler is written to files.
