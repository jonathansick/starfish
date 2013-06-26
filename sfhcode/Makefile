CC = gcc
FORT = f77
CFLAGS = -c
COMDIR = ../commoncode
BINDIR = ../.
HEADER = sfh.h
DOBJECTS = sfh.o amoeba.o dfstat.o direction.o errs.o \
fitstat.o grid.o logp.o simplex.o search.o step.o
OBJECTS = $(DOBJECTS) dchi.o random.o parseline.o utilities.o

sfh : $(OBJECTS) $(HEADER)
	$(FORT) -o sfh $(OBJECTS) 
	cp -f sfh $(BINDIR)

$(DOBJECTS) : $(HEADER)

random.o : $(COMDIR)/random.f
	$(FORT) $(CFLAGS) $(COMDIR)/random.f 
parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c 
utilities.o : $(COMDIR)/utilities.f
	$(FORT) $(CFLAGS) $(COMDIR)/utilities.f 

clean : 
	rm $(OBJECTS)
nuke :
	rm sfh ../sfh $(OBJECTS)

