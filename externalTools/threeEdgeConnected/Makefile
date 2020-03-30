rootPath = ../../
include ../../include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

all : ${LIBDIR}/3EdgeConnected.a ${BINDIR}/3EdgeTests

${LIBDIR}/3EdgeConnected.a : ${libSources} ${libHeaders} ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} -c ${libSources}
	${AR} rc 3EdgeConnected.a *.o
	${RANLIB} 3EdgeConnected.a 
	mv 3EdgeConnected.a ${LIBDIR}/
	cp ${libHeaders} ${LIBDIR}/

clean : 
	rm -f ${LIBDIR}/3EdgeConnected.a ${BINDIR}/3EdgeTests *.o

${BINDIR}/3EdgeTests : ${libTests} ${libSources} ${libHeaders} ${LIBDIR}/3EdgeConnected.a ${LIBDEPENDS}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/3EdgeTests ${libTests} ${LIBDIR}/3EdgeConnected.a ${LDLIBS}
