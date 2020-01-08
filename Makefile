#############################################################################
# Make file for AutoMix sampler
# Written by: David Hastie, University of Bristol
# 
# If you use the AutoMix sampler please reference the author and work
# as instructed on the website http://www.davidhastie.me.uk/AutoMix
# Initially this reference is likely to be the Ph.D. thesis that introduces 
# the AutoMix sampler. However, this will hopefully change to be a published
# paper in the not too distant future.
#
#############################################################################

# The file contains commands for compiling programs using optimisation and 
# no debugging.
# If debug information is needed (-g compiler flag) call `make DEBUG=1`
# on the command line.
# Flags and compiler names should be changed as necessary

ifdef DEBUG
CFLAGS=-g
else
CFLAGS=-O3
endif

# Libraries
LIBS=-lm

OBJS = automix.o logwrite.o utils.o

# Compiler flags:
#       -g      -- Enable debugging through gdb
#       -O3     -- Optimise progs 
#       -Wall   -- Turn on all warnings
#       -lm     -- Use the maths library
#       -c      -- Compile flag for dependencies

###### Type "make all" to make all files in the folder #####

all: amtoy1 amtoy2 amcpt amcptrs amrb9 amddi

###### Normal (already debugged) progs ############

# Toy example 1
am%: user%.o main.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ main.c $< $(OBJS) $(LIBS)

### AutoMix dependencies

# Rule for all object files
%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

### User supplied functions

# Rule for all user* object files
user%.o: user%.c user.h
	$(CC) $(CFLAGS) -c $<

# Rule for ddi (includes ddidata.h)
userddi.o: userddi.c user.h ddidata.h
	$(CC) $(CFLAGS) -c $<

###### Type "make clean" to remove all executables and object files ####

clean:
	- rm *.o
	- rm -r *.dSYM
	- rm am*
