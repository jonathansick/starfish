CC = gcc
FC = gfortran
CFLAGS = -c
FCFLAGS = -std=legacy
COMDIR = ../commoncode
BINDIR = ../.
GETEEP_OBJECTS = geteep.o 
INTERP_OBJECTS = interp.o 
COMMON_OBJECTS = parseline.o utilities.o

all     : geteep interp

geteep  : $(GETEEP_OBJECTS) $(COMMON_OBJECTS)
	$(FC) -o geteep $(GETEEP_OBJECTS) $(COMMON_OBJECTS)
	cp -f geteep $(BINDIR)

interp  : $(INTERP_OBJECTS) $(COMMON_OBJECTS)
	$(FC) -o interp $(INTERP_OBJECTS) $(COMMON_OBJECTS)
	cp -f interp $(BINDIR)

parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c

utilities.o : $(COMDIR)/utilities.f
	$(FC) $(FCFLAGS) -c $(COMDIR)/utilities.f

clean   :
	rm $(GETEEP_OBJECTS) $(INTERP_OBJECTS) $(COMMON_OBJECTS)

nuke    :
	rm geteep $(BINDIR)/geteep interp $(BINDIR)/interp $(GETEEP_OBJECTS) $(INTERP_OBJECTS) $(COMMON_OBJECTS)
