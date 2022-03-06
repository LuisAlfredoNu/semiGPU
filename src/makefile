.PHONY : all runQuickTest install update


SET_DEVICE :=$(DEVICE)

ifeq ($(strip $(DEVICE)),)
  SET_DEVICE := seq
endif

all:
	 @echo $(SET_DEVICE)
	@$(MAKE) -C common DEBUGVERSION=$(DEBUGVERSION) OPENACC=$(OPENACC) DEVICE=$(SET_DEVICE)
	@$(MAKE) -C test   DEBUGVERSION=$(DEBUGVERSION) OPENACC=$(OPENACC) DEVICE=$(SET_DEVICE)
	@$(MAKE) -C App    DEBUGVERSION=$(DEBUGVERSION) OPENACC=$(OPENACC) DEVICE=$(SET_DEVICE)

runQuickTest:
	@cd ../bin; ./alignandcomparemolecule.x ../bigTest/TesisExample/Opt_Geo_2/test-R_out_G09.xyz ../bigTest/TesisExample/Opt_Geo_2/test-S_out_G09.xyz

install:
	@cd ../bin; cp alignandcomparemolecule.x /usr/local/bin/alignandcomparemolecule
	@cd ../bin; cp rotationmolecules.x /usr/local/bin/rotationmolecules

clean:
	@cd common; make clean; 
	@cd test; make clean; 
	@cd App; make clean;

update: clean all
