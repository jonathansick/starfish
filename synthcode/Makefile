CC = gcc
FORT = f77
CFLAGS = -c
COMDIR = ../commoncode
BINDIR = ../.
OBJECTS = synth.o axes.o bin_dmags.o bin_interp.o bincombine.o \
crowd.o cumfcn.o delta_axes.o dmags.o fakecrowd.o mkcfn.o outpix.o \
pickphot.o readcrowd.o readphot.o redden.o random.o parseline.o \
utilities.o

synth : $(OBJECTS) synth.h
	$(FORT) -o synth $(OBJECTS) 
	cp -f synth $(BINDIR)

random.o : $(COMDIR)/random.f
	$(FORT) $(CFLAGS) $(COMDIR)/random.f
parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c 
utilities.o : $(COMDIR)/utilities.f
	$(FORT) $(CFLAGS) $(COMDIR)/utilities.f

clean : 
	rm $(OBJECTS)
nuke :
	rm synth ../synth $(OBJECTS)

