# Make file for MaxChi and Phylpro

CC = gcc
LIBRARIES = -L${R_HOME}/lib

OBJECT-FILES-MAXCHI = init.o getMaxChi.o MaxChi.o Rmaxchi.o
OBJECT-FILES-PHYLPRO = init.o getPhylpro.o Phylpro.o Rphylpro.o

all: MaxChi.so Phylpro.so

MaxChi.so: ${OBJECT-FILES-MAXCHI}
	${CC} -shared -o $@ ${OBJECT-FILES-MAXCHI} ${LIBRARIES}

Phylpro.so: ${OBJECT-FILES-PHYLPRO}
	${CC} -shared -o $@ ${OBJECT-FILES-PHYLPRO} ${LIBRARIES}

clean:
	rm init.o getMaxChi.o MaxChi.o getPhylpro.o Phylpro.o\
        Rphylpro.o Rmaxchi.o MaxChi.so Phylpro.so

