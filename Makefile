#############################################################################
# Make file for AutoMix sampler
# Written by: David Hastie, University of Bristol
# Last edited: 14/11/04
# 
# If you use the AutoMix sampler please reference the author and work
# as instructed on the website http://www.davidhastie.me.uk/AutoMix
# Initially this reference is likely to be the Ph.D. thesis that introduces 
# the AutoMix sampler. However, this will hopefully change to be a published
# paper in the not too distant future.
#
# All code has only been tested using the GNU compiler gcc.
#############################################################################

# The file contains commands for compiling programs using optimisation and 
# no debugging. Also included are commands for compiling with debugging.
# The executables and object files for the debugging are identified by 
# the final letter "d". For example, amtoy1.exe is the AutoMix sampler for the
# toy problem (see Ph.D. thesis section 5.5.1 for details) and amtoy1d.exe is 
# the version that allows debugging using a debugger such as gdb
# Flags and compiler names should be changed as necessary.

# Normal optimised flags
CFLAGS=-O3 -o

# Debugging flags
CFLAGSD=-g -o

# Dependency flags
DEPFLAGS=-c

# Dependency flags with debugging
DEPFLAGSD=-g -c

# Libraries
LIB=-lm

# Compiler flags:
#       -g      -- Enable debugging through gdb
#       -O3     -- Optimise progs 
#       -Wall   -- Turn on all warnings
#       -lm     -- Use the maths library
#       -c      -- Compile flag for dependencies


###### Type "make all" to make all files in the folder #####

all: amtoy1 amtoy2 amcpt amcptrs amrb9 amddi
alldebug: amtoy1d amtoy2d amcptd amcptrsd amrb9d amddid

###### Normal (already debugged) progs ############

# Toy example 1
amtoy1: main.c automix.o logwrite.o usertoy1.o utils.o
	$(CC) $(CFLAGS) amtoy1 main.c automix.o logwrite.o usertoy1.o utils.o $(LIB)

# Toy example 2
amtoy2: main.c automix.o logwrite.o usertoy2.o utils.o
	$(CC) $(CFLAGS) amtoy2 main.c automix.o logwrite.o usertoy2.o utils.o $(LIB)

# Change point problem
amcpt: main.c automix.o logwrite.o usercpt.o utils.o
	$(CC) $(CFLAGS) amcpt main.c automix.o logwrite.o usercpt.o utils.o $(LIB)

# Rescaled change point problem
amcptrs: main.c automix.o logwrite.o usercptrs.o utils.o
	$(CC) $(CFLAGS) amcptrs main.c automix.o logwrite.o usercptrs.o utils.o $(LIB)

# Rb9 problem
amrb9: main.c automix.o logwrite.o userrb9.o utils.o
	$(CC) $(CFLAGS) amrb9 main.c automix.o logwrite.o userrb9.o utils.o $(LIB)

# DDI Clinical trial problem
amddi: main.c automix.o logwrite.o ddidata.h userddi.o utils.o
	$(CC) $(CFLAGS) amddi main.c automix.o logwrite.o userddi.o utils.o $(LIB)

### AutoMix dependencies (already debugged)

# Utils
utils.o: utils.c utils.h
	$(CC) $(DEPFLAGS) utils.c -DDOUB -DRETS

# Log Write
logwrite.o: logwrite.c logwrite.h automix.h utils.h
	$(CC) $(DEPFLAGS) logwrite.c -DDOUB -DRETS

# AutoMix
automix.o: automix.c automix.h user.h
	$(CC) $(DEPFLAGS) automix.c -DDOUB -DRETS

### User supplied functions (already debugged)

# Toy example 1
usertoy1.o: usertoy1.c
	$(CC) $(DEPFLAGS) usertoy1.c

# Toy example 2
usertoy2.o: usertoy2.c
	$(CC) $(DEPFLAGS) usertoy2.c

# Change point problem
usercpt.o: usercpt.c
	$(CC) $(DEPFLAGS) usercpt.c

# Rescaled change point problem
usercptrs.o: usercptrs.c
	$(CC) $(DEPFLAGS) usercptrs.c

# Rb9 problem
userrb9.o: userrb9.c
	$(CC) $(DEPFLAGS) userrb9.c

# DDI clinical trial problem
userddi.o: userddi.c
	$(CC) $(DEPFLAGS) userddi.c

###### Progs to be debugged ############

# Toy example 1
amtoy1d: main.c automixd.o logwrited.o usertoy1d.o utilsd.o
	$(CC) $(CFLAGSD) amtoy1d main.c automixd.o logwrited.o usertoy1d.o utilsd.o $(LIB)

# Toy example 2
amtoy2d: main.c automixd.o logwrited.o usertoy2d.o utilsd.o
	$(CC) $(CFLAGSD) amtoy2d main.c automixd.o logwrited.o usertoy2d.o utilsd.o $(LIB)

# Change point problem
amcptd: main.c automixd.o logwrited.o usercptd.o utilsd.o
	$(CC) $(CFLAGSD) amcptd main.c automixd.o logwrited.o usercptd.o utilsd.o $(LIB)

# Rescaled change point problem
amcptrsd: main.c automixd.o logwrited.o usercptrsd.o utilsd.o
	$(CC) $(CFLAGSD) amcptrsd main.c automixd.o logwrited.o usercptrsd.o utilsd.o $(LIB)

# Rb9 problem
amrb9d: main.c automixd.o logwrited.o userrb9d.o utilsd.o
	$(CC) $(CFLAGSD) amrb9d main.c automixd.o logwrited.o userrb9d.o utilsd.o $(LIB)

# DDI Clinical trial problem
amddid: main.c automixd.o logwrited.o ddidata.h userddid.o utilsd.o
	$(CC) $(CFLAGSD) amddid main.c automixd.o logwrited.o userddid.o utilsd.o $(LIB)

# (Old) Toy problems compiled with automix2.c program, implementing
# adaptation through regeneration (automix2.c not included in distribution) 

### AutoMix dependencies (to be debugged)

# Utils
utilsd.o: utils.c utils.h
	$(CC) $(DEPFLAGSD) utils.c -DDOUB -DRETS -o utilsd.o

# Log Write
logwrited.o: logwrite.c logwrite.h automix.h utils.h
	$(CC) $(DEPFLAGSD) logwrite.c -DDOUB -DRETS -o logwrited.o

# AutoMix
automixd.o: automix.c automix.h user.h
	$(CC) $(DEPFLAGSD) automix.c -DDOUB -DRETS -o automixd.o


### User supplied functions (to be debugged)

# Toy example 1
usertoy1d.o: usertoy1.c
	$(CC) $(DEPFLAGSD) usertoy1.c -o usertoy1d.o

# Toy example 2
usertoy2d.o: usertoy2.c
	$(CC) $(DEPFLAGSD) usertoy2.c -o usertoy2d.o

# Change point problem
usercptd.o: usercpt.c
	$(CC) $(DEPFLAGSD) usercpt.c -o usercptd.o

# Rescaled change point problem
usercptrsd.o: usercptrs.c
	$(CC) $(DEPFLAGSD) usercptrs.c -o usercptrsd.o

# Rb9 problem
userrb9d.o: userrb9.c
	$(CC) $(DEPFLAGSD) userrb9.c -o userrb9d.o

# DDI clinical trial problem
userddid.o: userddi.c
	$(CC) $(DEPFLAGSD) userddi.c -o userddid.o

###### Type "make clean" to remove all executables and object files ####

clean:
	- rm *.o
	- rm -r *.dSYM
	- rm am*
