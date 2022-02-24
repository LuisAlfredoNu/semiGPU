.PHONY : all runQuickTest install update

all:
	@cd common; make -j 8; 
	@cd test; make -j 8;
	@cd App; make -j 8;

runQuickTest:
	@cd ../bin; ./alignandcomparemolecule.x ../bigTest/TesisExample/Opt_Geo_2/test-R_out_G09.xyz ../bigTest/TesisExample/Opt_Geo_2/test-S_out_G09.xyz

install:
	@cd ../bin; cp alignandcomparemolecule.x /usr/local/bin/alignandcomparemolecule
	@cd ../bin; cp rotationmolecules.x /usr/local/bin/rotationmolecules

clean:
	@cd common; make clean; 
	@cd test; make clean; 

update: clean all
