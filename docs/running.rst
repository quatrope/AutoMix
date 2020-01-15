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

  - **-m D** controls the mode of the sampler.
  `D=0` is mixture fitting;
  `D=1` skips stage 1 and 2 if a file containing the mixture parameters is supplied;
  `D=2` fits AutoMix version of AutoRJ sampler (see Green, 2003 - full reference in thesis).
  Default uses `D=0`.
  - **-n I** run the sampler for `max(I,nkk*10000,100000)` iterations in the stage 1 RWM for each model k. (Default uses I=100000)
  - **-N I** run the sampler for I Reversible jump iterations in stage 3. (Default uses I=100000).
  - **-s I** initialises the random number generator with seed I. (Default uses clock as seed).
  - **-a A** controls whether or not adaptation is done in stage 3 RJ. If `A=0` no adaptation is done, if `A=1` adaptation is done. (Default has `A=1`).
  - **-p P** controls whether or not random permutation is done in stage 3 RJ. If `P=0` no permutation is done, if `P=1` permutation is done. Default has `P=0`.
  - **-t I**      Controls whether standard Normal or t distributed variables are used in RWM and in RJ moves. If `I=0` Normal variables are used, otherwise t-distributed variables with I degrees of freedom are used. Default `I=0`.
  - **-f F** Uses the string F as the bases for filenames (e.g. if `F=output`, filenames are `output_log.data`, `output_mix.data` etc). (Default is `F=output`)
  - **-h, --help** Prints help information on command line arguments and exits.
  - **-v, --version** Prints version number and exits.

Example
^^^^^^^

As an example, typing::

    amtoy1 -m 0 -N 1000000 -p 1 -f toy1

runs the optimized mixture fitting version of the toy1 problem (see thesis, section 5.5.1) with 1 million RJ sweeps, enabling permutation and storing the output in files of the type `toy1_***.data`.
Running the sampler produces a summary of how the run is progressing.

For each of the models:

* In **stage 1** a countdown of the number of iterations remaining is printed to screen;
* In **stage 2** a summary of the mixture fitting is printed to screen.
  This summary consists of a countdown of the number of components in the current mixture, with the iteration number that the last component 
  was removed and an indicator `n` if the component was annihilated naturally, and `f` if the annihilation was forced.
* In the RJ **stage 3** a countdown of the number of iterations remaining is printed to screen. 
  No summary statistics are printed to screen. Instead all output from the sampler is written to files.
