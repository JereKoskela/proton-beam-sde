# For iOS only
GSL_PREFIX = $(shell brew --prefix gsl)
LCONFIG_PREFIX = $(shell brew --prefix libconfig)

CFLAGS=-Wall -Wextra -O3 -I$(GSL_PREFIX)/include -I$(LCONFIG_PREFIX)/include
LDFLAGS= -L$(GSL_PREFIX)/lib -L$(LCONFIG_PREFIX)/lib -lgsl -lgslcblas -lconfig++
CC = g++

main: simulate.cc
	${CC} ${CFLAGS} -o simulate simulate.cc ${LDFLAGS}
