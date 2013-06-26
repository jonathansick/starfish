CC = gcc
F77 = f77
CFLAGS = -c
COMDIR = ../commoncode
BINDIR = ../.
GETEEP_OBJECTS = geteep.o 
INTERP_OBJECTS = interp.o 
COMMON_OBJECTS = parseline.o utilities.o

all     : geteep interp

geteep  : $(GETEEP_OBJECTS) $(COMMON_OBJECTS)
	$(F77) -o geteep $(GETEEP_OBJECTS) $(COMMON_OBJECTS)
	cp -f geteep $(BINDIR)

interp  : $(INTERP_OBJECTS) $(COMMON_OBJECTS)
	$(F77) -o interp $(INTERP_OBJECTS) $(COMMON_OBJECTS)
	cp -f interp $(BINDIR)

parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c

utilities.o : $(COMDIR)/utilities.f
	$(F77) $(CFLAGS) $(COMDIR)/utilities.f

clean   :
	rm $(GETEEP_OBJECTS) $(INTERP_OBJECTS) $(COMMON_OBJECTS)

nuke    :
	rm geteep $(BINDIR)/geteep interp $(BINDIR)/interp $(GETEEP_OBJECTS) $(INTERP_OBJECTS) $(COMMON_OBJECTS)
