CC = gcc
FC = gfortran
CFLAGS = -c
FCFLAGS = -std=legacy
COMDIR = ../commoncode
BINDIR = ../.
HEADER = sfh.h
DOBJECTS = sfh.o amoeba.o dfstat.o direction.o errs.o \
fitstat.o grid.o logp.o simplex.o search.o step.o
OBJECTS = $(DOBJECTS) dchi.o random.o parseline.o message.o utilities.o

sfh : $(OBJECTS) $(HEADER)
	$(FC) -o sfh $(OBJECTS) 
	cp -f sfh $(BINDIR)

sfh.o :
	$(FC) $(FCFLAGS) -c sfh.f
amoeba.o :
	$(FC) $(FCFLAGS) -c amoeba.f
dfstat.o :
	$(FC) $(FCFLAGS) -c dfstat.f
direction.o :
	$(FC) $(FCFLAGS) -c direction.f
errs.o :
	$(FC) $(FCFLAGS) -c errs.f
fitstat.o :
	$(FC) $(FCFLAGS) -c fitstat.f
grid.o :
	$(FC) $(FCFLAGS) -c grid.f
simplex.o :
	$(FC) $(FCFLAGS) -c simplex.f
search.o :
	$(FC) $(FCFLAGS) -c search.f
step.o :
	$(FC) $(FCFLAGS) -c step.f
dchi.o :
	$(FC) $(FCFLAGS) -c dchi.f

random.o : $(COMDIR)/random.f
	$(FC) $(FCFLAGS) -c $(COMDIR)/random.f 
parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c 
message.o : $(COMDIR)/message.c
	$(CC) $(CFLAGS) $(COMDIR)/message.c 
utilities.o : $(COMDIR)/utilities.f
	$(FC) $(FCFLAGS) -c $(COMDIR)/utilities.f 

clean : 
	rm $(OBJECTS)
nuke :
	rm sfh ../sfh $(OBJECTS)

