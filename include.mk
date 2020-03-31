#Modify this variable to set the location of sonLib
sonLibRootDir ?= ${rootPath}/../sonLib
sonLibDir=${sonLibRootDir}/lib
#Use sonLib bin and lib dirs
BINDIR ?= ${sonLibRootDir}/bin
LIBDIR ?= ${sonLibDir}

include  ${sonLibRootDir}/include.mk

CPPFLAGS += -I${sonLibRootDir}/lib
LDLIBS = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a ${dblibs}
LIBDEPENDS = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a 
