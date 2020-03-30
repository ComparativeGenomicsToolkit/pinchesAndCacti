rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c
testBin = tests/testBin

all: all_libs all_progs
all_libs: externalToolsM ${LIBDIR}/stPinchesAndCacti.a
all_progs: all_libs
	${MAKE} ${BINDIR}/stPinchesAndCactiTests

externalToolsM : 
	cd externalTools && ${MAKE} all

${LIBDIR}/stPinchesAndCacti.a:  ${libSources} ${libHeaders} ${LIBDEPENDS} externalToolsM
	${CC} ${CPPFLAGS} ${CFLAGS} -c ${libSources}
	${AR} rc stPinchesAndCacti.a *.o
	${RANLIB} stPinchesAndCacti.a 
	mv stPinchesAndCacti.a ${LIBDIR}/
	cp ${libHeaders} ${LIBDIR}/

${BINDIR}/stPinchesAndCactiTests : ${libTests} ${libSources} ${libHeaders} ${LIBDEPENDS} externalToolsM ${LIBDIR}/3EdgeConnected.a
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stPinchesAndCactiTests ${libTests} ${libSources} ${LIBDIR}/3EdgeConnected.a ${LDLIBS}

clean : 
	cd externalTools && ${MAKE} clean
	rm -f *.o
	rm -f ${LIBDIR}/stPinchesAndCacti.a ${BINDIR}/stPinchesAndCactiTests

test : all
	${BINDIR}/stPinchesAndCactiTests
