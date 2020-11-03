.PHONY : all runQuickTest install update

all:
	@cd alignandcomparemolecule; make; 
	@cd rotationmolecules; make;
	@cd ..; mkdir -p bin; cp src/alignandcomparemolecule/alignandcomparemolecule.x bin/
	@cd ..; mkdir -p bin; cp src/rotationmolecules/rotationmolecules.x bin/

runQuickTest:
	@cd ../bin; ./alignandcomparemolecule.x ../bigTest/TesisExample/Opt_Geo_2/test-R_out_G09.xyz ../bigTest/TesisExample/Opt_Geo_2/test-S_out_G09.xyz

install:
	@cd ../bin; cp alignandcomparemolecule.x /usr/local/bin/alignandcomparemolecule
	@cd ../bin; cp rotationmolecules.x /usr/local/bin/rotationmolecules

clean:
	@cd alignandcomparemolecule; make fullclean; 

	@cd alignandcomparemolecule; make fullclean; 
update: clean all
