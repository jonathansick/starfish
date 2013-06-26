FORT = f77
CC = gcc
CFLAGS = -c
SYNTHDIR = ../synthcode
COMDIR = ../commoncode
OBJECTS = repop.o axes.o bin_dmags.o bin_interp.o bincombine.o \
crowd.o cumfcn.o delta_axes.o dmags.o fakecrowd.o mkcfn.o outpix.o \
pickphot.o readcrowd.o readphot.o redden.o random.o parseline.o \
utilities.o

repop : $(OBJECTS) $(SYNTHDIR)/synth.h
	$(FORT) -o repop $(OBJECTS) 

repop.o : repop.f
	$(FORT) $(CFLAGS) repop.f

axes.o : $(SYNTHDIR)/axes.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/axes.f

bin_dmags.o : $(SYNTHDIR)/bin_dmags.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/bin_dmags.f

bin_interp.o : $(SYNTHDIR)/bin_interp.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/bin_interp.f

bincombine.o : $(SYNTHDIR)/bincombine.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/bincombine.f

crowd.o : $(SYNTHDIR)/crowd.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/crowd.f

cumfcn.o : $(SYNTHDIR)/cumfcn.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/cumfcn.f

delta_axes.o : $(SYNTHDIR)/delta_axes.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/delta_axes.f

dmags.o : $(SYNTHDIR)/dmags.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/dmags.f

fakecrowd.o : $(SYNTHDIR)/fakecrowd.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/fakecrowd.f

mkcfn.o : $(SYNTHDIR)/mkcfn.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/mkcfn.f

outpix.o : $(SYNTHDIR)/outpix.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/outpix.f

pickphot.o : $(SYNTHDIR)/pickphot.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/pickphot.f

readcrowd.o : $(SYNTHDIR)/readcrowd.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/readcrowd.f

readphot.o : $(SYNTHDIR)/readphot.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/readphot.f

redden.o : $(SYNTHDIR)/redden.f
	$(FORT) $(CFLAGS) $(SYNTHDIR)/redden.f

random.o : $(COMDIR)/random.f
	$(FORT) $(CFLAGS) $(COMDIR)/random.f

parseline.o : $(COMDIR)/parseline.c
	$(CC) $(CFLAGS) $(COMDIR)/parseline.c 

utilities.o : $(COMDIR)/utilities.f
	$(FORT) $(CFLAGS) $(COMDIR)/utilities.f

clean : 
	rm $(OBJECTS)
nuke :
	rm repop $(OBJECTS)

