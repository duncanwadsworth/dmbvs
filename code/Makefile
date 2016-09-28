## Makefile for dmbvs

## file and path names
## LIBS and INCS may need to be changed
INPUT = dmbvs
LIBS = /opt/local/lib
INCS = /opt/local/include

compile:
	#clang -I$(INCS) -L$(LIBS) -Wall $(INPUT).c -o $(INPUT).x -lgsl -lgslcblas -lm
	#/local/gcc-5.2.0/bin/gcc -I/local/gsl-2.1/include -L/local/gsl-2.1/lib -Wall $(INPUT).c -o $(INPUT).x -lgsl -lgslcblas -lm
	gcc -I$(INCS) -L$(LIBS) -Wall $(INPUT).c -o $(INPUT).x -lgsl -lgslcblas -lm

clean:
	rm -f $(INPUT).x
	if test -a $(INPUT).o ; then rm -f $(INPUT).o ; fi

