.PHONY : all runQuickTest install update

SET_DEVICE := $(DEVICE)

ifeq ($(strip $(DEVICE)),)
  SET_DEVICE := seq
endif

all:
	 @echo $(SET_DEVICE)
	@$(MAKE) -C common DEBUGVERSION=$(DEBUGVERSION) OPENACC=$(OPENACC) DEVICE=$(SET_DEVICE)
	@$(MAKE) -C test   DEBUGVERSION=$(DEBUGVERSION) OPENACC=$(OPENACC) DEVICE=$(SET_DEVICE)
	@$(MAKE) -C App    DEBUGVERSION=$(DEBUGVERSION) OPENACC=$(OPENACC) DEVICE=$(SET_DEVICE)

runQuickTest:
	@cd test; bash fullTester.bsh ../Geometries4Test/03_formaldehyde.xyz

runParallelTest:
	@cd test; bash fullTester.bsh ../Geometries4Test/03_formaldehyde.xyz 1

install:
	@cd ../bin; cp alignandcomparemolecule.x /usr/local/bin/alignandcomparemolecule
	@cd ../bin; cp rotationmolecules.x /usr/local/bin/rotationmolecules

clean:
	@cd common; make clean; 
	@cd test; make clean; 
	@cd App; make clean;

fullclean:
	@cd common; make clean; 
	@cd test; make fullclean; 
	@cd App; make fullclean;


update: clean all
