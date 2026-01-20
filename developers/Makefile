#
# Tests the installed version of meddly
# (headers and compiled library)
#

all: nqueens

clean:
	rm nqueens nqueens.o

nqueens: nqueens.o
	$(CXX) -L../lib nqueens.o -lmeddly -o nqueens

nqueens.o: nqueens.cc
	$(CXX) -c -I../include/meddly nqueens.cc
