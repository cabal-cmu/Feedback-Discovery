include ${FSLCONFDIR}/default.mk

PROJNAME = nma

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} -I${INC_ZLIB} -I${INC_BOOST}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_ZLIB}

LIBS =  -lutils -lmiscplot -lmiscpic -lwarpfns -lnewimage -lmiscmaths -lm -lnewmat -lfslio -lniftiio -lznz -lprob -lgdc -lgd -lpng -lz 

RUNTCLS = Nma

SCRIPTS = nma_launch

XFILES = nma nma_generate

OBJS = nma_generate_options.o nmaoptions.o nma_manager.o model.o subject.o node.o halfcos_hrf.o mcmc_mh.o mcmc_component.o vb_component.o 

#OPTFLAGS =

all:	${XFILES}

nma:${OBJS} nma.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} nma.o ${LIBS}

nma_generate:${OBJS} nma_generate.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} nma_generate.o ${LIBS}
