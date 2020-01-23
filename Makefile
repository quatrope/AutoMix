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
CFLAGS=-O3 -Wall
endif

ifdef PREFIX
LIB_INST_DIR=$(PREFIX)/lib
INC_INST_DIR=$(PREFIX)/include
else
LIB_INST_DIR=/usr/local/lib
INC_INST_DIR=/usr/local/include
endif

# Libraries
LIBS=-lm
LIBOBJS=automix.o
SRC_DIR=src
LIB_DIR=$(SRC_DIR)/libautomix
EXMP_DIR=$(SRC_DIR)/user_examples
INC_DIRS=$(EXMP_DIR) $(LIB_DIR)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

all: libautomix.so

examples: amtoy1 amtoy2 amcpt amcptrs amrb9 amddi

libautomix.so: $(LIBOBJS)
	$(CC) -shared -o libautomix.so $(LIBOBJS) $(LIBS)

###### EXAMPLE PROGRAMS ############

# Toy example 1
am%: user%.o logwrite.o $(EXMP_DIR)/main.c libautomix.so
	$(CC) $(CFLAGS) -o $@ $< $(EXMP_DIR)/main.c logwrite.o -L./ -lautomix $(LIBS) $(INC_FLAGS)

### AutoMix dependencies

# Rule for all object files
%.o: $(LIB_DIR)/%.c $(LIB_DIR)/%.h
	$(CC) $(CFLAGS) -fPIC -c $<

### User supplied functions

# Rule for all user* object files
user%.o: $(EXMP_DIR)/user%.c $(EXMP_DIR)/user.h
	$(CC) $(CFLAGS) -c $<

# Rule for ddi (includes ddidata.h)
userddi.o: $(EXMP_DIR)/userddi.c $(EXMP_DIR)/user.h $(EXMP_DIR)/ddidata.h
	$(CC) $(CFLAGS) -c $<

logwrite.o: $(EXMP_DIR)/logwrite.c $(EXMP_DIR)/logwrite.h
	$(CC) $(CFLAGS) -c $< -I$(LIB_DIR)

tutorial: $(EXMP_DIR)/tutorial.c libautomix.so
	$(CC) $(CFLAGS) -o $@ $< -L./ -lautomix $(LIBS) $(INC_FLAGS)

###### Type "make clean" to remove all executables and object files ####

test: tests/main.c logwrite.o libautomix.so
	$(CC) $< logwrite.o -L./ -lautomix -lm -I$(LIB_DIR) -I$(EXMP_DIR) -o $@

install:
	- mkdir $(LIB_INST_DIR)/
	- mkdir $(INC_INST_DIR)/
	cp libautomix.so $(LIB_INST_DIR)/
	cp $(LIB_DIR)/automix.h $(INC_INST_DIR)/

clean:
	- rm *.o
	- rm amtoy1 amtoy2 amcpt amcptrs amrb9 amddi
	- rm test
	- rm tutorial
