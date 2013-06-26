FORT = f77
CC = gcc
CFLAGS = -c
SFHDIR = ../sfhcode
COMDIR = ../commoncode
BINDIR = ../.
OBJECTS = testchi.o fitstat.o grid.o parseline.o utilities.o 

testchi     : $(OBJECTS) $(SFHDIR)/sfh.h
	$(FORT) -o testchi $(OBJECTS) 
	cp -f testchi $(BINDIR)/testchi

fitstat.o   : $(SFHDIR)/fitstat.f
	$(FORT) $(CFLAGS) $(SFHDIR)/fitstat.f

grid.o      : $(SFHDIR)/grid.f
	$(FORT) $(CFLAGS) $(SFHDIR)/grid.f

parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c 

utilities.o : $(COMDIR)/utilities.f
	$(FORT) $(CFLAGS) $(COMDIR)/utilities.f

clean : 
	rm $(OBJECTS)
nuke :
	rm testchi $(BINDIR)/testchi $(OBJECTS)
