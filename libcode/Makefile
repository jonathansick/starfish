# The compiler
FC = gfortran
CC = gcc
CFLAGS = -c
FCFLAGS = -std=legacy -O3
COMDIR = ../commoncode
BINDIR = ../.
HEADER = lib.h
DOBJECTS = lib.o
OBJECTS = $(DOBJECTS) polfit.o pvalue.o parseline.o utilities.o

mklib       : $(OBJECTS) $(HEADER)
	$(FC) -o mklib $(OBJECTS)
	cp -f mklib $(BINDIR)

lib.o :
	$(FC) $(FCFLAGS) -c lib.f
polfit.o :
	$(FC) $(FCFLAGS) -c polfit.f
pvalue.o :
	$(FC) $(FCFLAGS) -c pvalue.f

parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c

utilities.o : $(COMDIR)/utilities.f
	$(FC) $(FCFLAGS) -c $(COMDIR)/utilities.f

clean       : 
		rm $(OBJECTS)
nuke        :
		rm mklib $(BINDIR)/mklib $(OBJECTS)
