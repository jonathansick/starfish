CC = gcc
FC = gfortran
CFLAGS = -c
FCFLAGS = -std=legacy -O3
COMDIR = ../commoncode
BINDIR = ../.
HEADER = synth.h
DOBJECTS = synth.o mkcfn.o axes.o bin_dmags.o bin_interp.o bincombine.o \
crowd.o cumfcn.o delta_axes.o dmags.o fakecrowd.o outpix.o \
pickphot.o readcrowd.o readphot.o redden.o scatter.o 
OBJECTS = $(DOBJECTS) random.o parseline.o utilities.o

synth : $(OBJECTS) $(HEADER)
	$(FC) $(FCFLAGS) -o synth $(OBJECTS) 
	cp -f synth $(BINDIR)

synth.o : 
	$(FC) $(FCFLAGS) -c synth.f
axes.o : 
	$(FC) $(FCFLAGS) -c axes.f
bin_dmags.o : 
	$(FC) $(FCFLAGS) -c bin_dmags.f
bin_interp.o : 
	$(FC) $(FCFLAGS) -c bin_interp.f
bincombine.o : 
	$(FC) $(FCFLAGS) -c bincombine.f
crowd.o : 
	$(FC) $(FCFLAGS) -c crowd.f
cumfcn.o : 
	$(FC) $(FCFLAGS) -c cumfcn.f
delta_axes.o : 
	$(FC) $(FCFLAGS) -c delta_axes.f
dmags.o : 
	$(FC) $(FCFLAGS) -c dmags.f
fakecrowd.o : 
	$(FC) $(FCFLAGS) -c fakecrowd.f
mkfcn.o : 
	$(FC) $(FCFLAGS) -c mkfcn.f
outpix.o : 
	$(FC) $(FCFLAGS) -c outpix.f
pickphot.o : 
	$(FC) $(FCFLAGS) -c pickphot.f
readcrowd.o : 
	$(FC) $(FCFLAGS) -c readcrowd.f
readphot.o : 
	$(FC) $(FCFLAGS) -c readphot.f
redden.o : 
	$(FC) $(FCFLAGS) -c redden.f
scatter.o : 
	$(FC) $(FCFLAGS) -c scatter.f

random.o : $(COMDIR)/random.f
	$(FC) $(FCFLAGS) $(CFLAGS) $(COMDIR)/random.f
parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c 
utilities.o : $(COMDIR)/utilities.f
	$(FC) $(FCFLAGS) -c $(COMDIR)/utilities.f

clean : 
	rm $(OBJECTS)
nuke :
	rm synth ../synth $(OBJECTS)

