COMPILER =g++
IDIR =./include
CXXFLAGS =-std=c++11 -Wall -O2 -I$(IDIR)
EXEC =full-model
LINK = $(CXX) $(CXXFLAGS)
COMPILE = $(CXX) $(LIBS) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -c

all: full-model snap2ascii summ2snap

full-model: full-model.o
	$(LINK) -o full-model full-model.o
snap2ascii: snap2ascii.o
	$(LINK) -o snap2ascii snap2ascii.o
summ2snap: summ2snap.o
	$(LINK) -o summ2snap summ2snap.o

full-model.o: ./src/evolve.cpp ./src/Gene.h ./src/Cell.h ./src/global.h
	$(COMPILE) -o full-model.o ./src/evolve.cpp
snap2ascii.o: ./tools/snap2ascii.cpp ./src/global.h
	$(COMPILE) -o snap2ascii.o ./tools/snap2ascii.cpp
summ2snap.o: ./tools/summ2snap.cpp ./src/PolyCell.h 
	$(COMPILE) -o summ2snap.o ./tools/summ2snap.cpp

clean:
	rm -f *.o