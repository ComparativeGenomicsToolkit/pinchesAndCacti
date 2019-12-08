rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c
testBin = tests/testBin

all: all_libs all_progs
all_libs: externalToolsM ${libPath}/stPinchesAndCacti.a
all_progs: all_libs
	${MAKE} ${binPath}/stPinchesAndCactiTests

externalToolsM : 
	cd externalTools && ${MAKE} all

${libPath}/stPinchesAndCacti.a:  ${libSources} ${libHeaders} ${basicLibsDependencies} externalToolsM
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc stPinchesAndCacti.a *.o
	ranlib stPinchesAndCacti.a 
	rm *.o
	mv stPinchesAndCacti.a ${libPath}/
	cp ${libHeaders} ${libPath}/

${binPath}/stPinchesAndCactiTests : ${libTests} ${libSources} ${libHeaders} ${basicLibsDependencies} externalToolsM ${libPath}/3EdgeConnected.a
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/stPinchesAndCactiTests ${libTests} ${libSources} ${basicLibs}  ${libPath}/3EdgeConnected.a

clean : 
	cd externalTools && ${MAKE} clean
	rm -f *.o
	rm -f ${libPath}/stPinchesAndCacti.a ${binPath}/stPinchesAndCactiTests

test : all
	${binPath}/stPinchesAndCactiTests
