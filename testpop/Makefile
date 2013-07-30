CC = gcc
FC = gfortran
CFLAGS = -c
FCFLAGS = -std=legacy
SYNTHDIR = ../synthcode
COMDIR = ../commoncode
OBJECTS = testpop.o axes.o bin_dmags.o bin_interp.o bincombine.o \
crowd.o cumfcn.o delta_axes.o dmags.o fakecrowd.o mkcfn.o outpix.o \
pickphot.o readcrowd.o readphot.o redden.o random.o parseline.o \
utilities.o

testpop : $(OBJECTS) $(SYNTHDIR)/synth.h
	$(FC) -o testpop $(OBJECTS) 

testpop.o : testpop.f
	$(FC) $(FCFLAGS) -c testpop.f

axes.o : $(SYNTHDIR)/axes.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/axes.f

bin_dmags.o : $(SYNTHDIR)/bin_dmags.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/bin_dmags.f

bin_interp.o : $(SYNTHDIR)/bin_interp.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/bin_interp.f

bincombine.o : $(SYNTHDIR)/bincombine.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/bincombine.f

crowd.o : $(SYNTHDIR)/crowd.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/crowd.f

cumfcn.o : $(SYNTHDIR)/cumfcn.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/cumfcn.f

delta_axes.o : $(SYNTHDIR)/delta_axes.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/delta_axes.f

dmags.o : $(SYNTHDIR)/dmags.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/dmags.f

fakecrowd.o : $(SYNTHDIR)/fakecrowd.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/fakecrowd.f

mkcfn.o : $(SYNTHDIR)/mkcfn.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/mkcfn.f

outpix.o : $(SYNTHDIR)/outpix.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/outpix.f

pickphot.o : $(SYNTHDIR)/pickphot.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/pickphot.f

readcrowd.o : $(SYNTHDIR)/readcrowd.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/readcrowd.f

readphot.o : $(SYNTHDIR)/readphot.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/readphot.f

redden.o : $(SYNTHDIR)/redden.f
	$(FC) $(FCFLAGS) -c $(SYNTHDIR)/redden.f

random.o : $(COMDIR)/random.f
	$(FC) $(FCFLAGS) -c $(COMDIR)/random.f

parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c 

utilities.o : $(COMDIR)/utilities.f
	$(FC) $(FCFLAGS) -c $(COMDIR)/utilities.f

clean : 
	rm $(OBJECTS)
nuke :
	rm testpop $(OBJECTS)

