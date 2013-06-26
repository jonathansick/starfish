CC = gcc
FORT = f77
CFLAGS = -c
COMDIR = ../commoncode
BINDIR = ../.
OBJECTS = sfh.o amoeba.o dchi.o dfstat.o direction.o \
errs.o fitstat.o grid.o logp.o simplex.o search.o step.o \
random.o parseline.o utilities.o #\
#ubbgrid.o bvvgrid.o viigrid.o 
# To use a custom CMD grid, replace the above files with 
# your custom grid subroutines (see grid.f).

sfh : $(OBJECTS) sfh.h
	$(FORT) -o sfh $(OBJECTS) 
	cp -f sfh $(BINDIR)

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

