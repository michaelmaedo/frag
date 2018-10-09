#========================================================================#
#                                                                        #
#  Makefile:    G M S H 2 C B                                            #
#               Mark_I                                                   #
#                                                                        #
#========================================================================#

# System Configuration Section
# ============================

# ** Shell **
SHELL = /bin/sh

# ** Suffixes **
.SUFFIXES = .f90 .o .mod

# ** Specifies a list of directories that 'make' should search **
#VPATH = src:src/utilities

# ** Compiler **
FC = gfortran -cpp -finit-local-zero -fmax-stack-var-size=30000 -O3\
	-std=f2003 -Wall -Wno-conversion

# ** Libraries **
#LIBS = -lblas

# ** Header Files **
INCLUDES = -I/usr/include

# ** Debug **
CDEBUG = -g


# ** Flags **
CFLAGS = -fdefault-double-8 -fno-range-check -Wall

# End of System Configuration Section
# ===================================
OBJS =  precision_mod.o parameters_mod.o jacob_mod.o mesh_frag_mod.o reorder_mod.o main.o

# ** Rules **
.PHONY: all
all: start

start:
	@echo "\n"
	@echo  "     F R A G"
	@echo  "     version 0"
	@echo  ""
	@echo  "     Insert finite elements with high aspect ratio,"
	@echo  "     also  called  interface  elements, between the"
	@echo  "     standard (bulk) elements of the original mesh"
	@echo  "     "
	@echo  ""
	@echo  "     Consulting : mmaedo@tamu.edu"
	@echo  ""
	@echo  ""
	@$(MAKE) frag
	@$(MAKE) clean
	@echo "\n * * * FRAG was successfully built * * *\n"

frag: $(OBJS)
	$(FC) $(CFLAGS) -o $@ $^

$(OBJS): $(FSRC)

FSRC:
	@FSRC = $(FSRC)

%.o %.mod: %.f90
	$(FC) -g -c $^

.PHONY: clean
clean:
	rm -f *.o *.mod

run:
	./gmsh2cb

# End of the makefile
