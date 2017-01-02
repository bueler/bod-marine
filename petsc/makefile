ALL: marine

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

bodvardsson : bodvardsson.o exactsolns.o chkopts
	-${CLINKER} -o bodvardsson bodvardsson.o exactsolns.o ${PETSC_LIB}
	${RM} bodvardsson.o exactsolns.o

marine : marine.o exactsolns.o chkopts
	-${CLINKER} -o marine marine.o exactsolns.o ${PETSC_LIB}
	${RM} marine.o exactsolns.o

