ALL: marine

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

exactTest.o : exactTestN.c

bodvardsson : bodvardsson.o exactTestN.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB} exactTestN.o

marine : marine.o exactTestN.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB} exactTestN.o

ex5: ex5.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}

ex9: ex9.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}

