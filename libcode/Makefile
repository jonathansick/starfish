# The compiler
FC = gfortran
CC = gcc
CFLAGS = -c
COMDIR = ../commoncode
BINDIR = ../.
HEADER = lib.h
DOBJECTS = lib.o
OBJECTS = $(DOBJECTS) polfit.o pvalue.o parseline.o utilities.o

mklib       : $(OBJECTS) $(HEADER)
	$(FC) -o mklib $(OBJECTS)
	cp -f mklib $(BINDIR)

$(DOBJECTS) : $(HEADER)

parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c

utilities.o : $(COMDIR)/utilities.f
	$(FC) $(CFLAGS) $(COMDIR)/utilities.f

clean       : 
		rm $(OBJECTS)
nuke        :
		rm mklib $(BINDIR)/mklib $(OBJECTS)
