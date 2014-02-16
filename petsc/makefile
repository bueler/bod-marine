ALL: marine

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

exactsolns.o : exactsolns.c

bodvardsson : bodvardsson.o exactTestN.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB} exactsolns.o

marine : marine.o exactsolns.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB} exactsolns.o

ex5: ex5.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}

ex9: ex9.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB}

