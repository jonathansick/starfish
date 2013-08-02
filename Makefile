# Global makefile for StarFISH

.PHONY: all

all: mklib synth sfh interp testchi repop/repop testpop/testpop

mklib:
	$(MAKE) --directory=libcode
synth:
	$(MAKE) --directory=synthcode
sfh:
	$(MAKE) --directory=sfhcode
interp:
	$(MAKE) --directory=interpcode
testchi:
	$(MAKE) --directory=testchicode
repop/repop:
	$(MAKE) --directory=repop
testpop/testpop:
	$(MAKE) --directory=testpop
