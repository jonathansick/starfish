CC = gcc
F77 = f77
CFLAGS = -c
COMDIR = ../commoncode
BINDIR = ../.
OBJECTS = lib.o polfit.o pvalue.o parseline.o utilities.o

mklib       : $(OBJECTS) lib.h
	$(F77) -o mklib $(OBJECTS)
	cp -f mklib $(BINDIR)

parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c

utilities.o : $(COMDIR)/utilities.f
	$(F77) $(CFLAGS) $(COMDIR)/utilities.f

clean       : 
		rm $(OBJECTS)
nuke        :
		rm mklib $(BINDIR)/mklib $(OBJECTS)
