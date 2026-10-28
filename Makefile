rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c
testBin = tests/testBin

all: all_libs all_progs
all_libs: ${LIBDIR}/stPinchesAndCacti.a
all_progs: all_libs
	${MAKE} ${BINDIR}/stPinchesAndCactiTests

${LIBDIR}/stPinchesAndCacti.a:  ${libSources} ${libHeaders} ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} -c ${libSources}
	${AR} rc stPinchesAndCacti.a *.o
	${RANLIB} stPinchesAndCacti.a 
	mv stPinchesAndCacti.a ${LIBDIR}/
	cp ${libHeaders} ${LIBDIR}/

${BINDIR}/stPinchesAndCactiTests : ${libTests} ${libSources} ${libHeaders} ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/stPinchesAndCactiTests ${libTests} ${libSources} ${LDLIBS}

clean : 
	rm -f *.o
	rm -f ${LIBDIR}/stPinchesAndCacti.a ${BINDIR}/stPinchesAndCactiTests

test : all
	${BINDIR}/stPinchesAndCactiTests
