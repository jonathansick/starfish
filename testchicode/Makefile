CC = gcc
FC = gfortran
CFLAGS = -c
FCFLAGS = -std=legacy
SFHDIR = ../sfhcode
COMDIR = ../commoncode
BINDIR = ../.
OBJECTS = testchi.o fitstat.o grid.o parseline.o utilities.o 

testchi     : $(OBJECTS) $(SFHDIR)/sfh.h
	$(FC) -o testchi $(OBJECTS) 
	cp -f testchi $(BINDIR)/testchi

fitstat.o   : $(SFHDIR)/fitstat.f
	$(FC) $(FCFLAGS) -c $(SFHDIR)/fitstat.f

grid.o      : $(SFHDIR)/grid.f
	$(FC) $(FCFLAGS) -c $(SFHDIR)/grid.f

parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c 

utilities.o : $(COMDIR)/utilities.f
	$(FC) $(FCFLAGS) -c $(COMDIR)/utilities.f

clean : 
	rm $(OBJECTS)
nuke :
	rm testchi $(BINDIR)/testchi $(OBJECTS)
