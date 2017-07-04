RM =rm
CP =cp
CC =g++

IDIR =./include
CXXFLAGS =-std=c++11 -Wall -O2 -I$(IDIR)
EXEC =sodapop
LINK = $(CXX) $(CXXFLAGS)
COMPILE = $(CXX) $(LIBS) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -c

SODAPOP = sodapop
SNAP2ASCII = snap2ascii
SUMM2SNAP = summ2snap



INSTALLDIR = /usr/local/bin

all: $(SODAPOP) $(SNAP2ASCII) $(SUMM2SNAP)
install:
	$(CP) $(SODAPOP) $(INSTALLDIR)/
	$(CP) $(SNAP2ASCII) $(INSTALLDIR)/
	$(CP) $(SUMM2SNAP) $(INSTALLDIR)/

$(SODAPOP): sodapop.o
	$(LINK) -o sodapop sodapop.o
$(SNAP2ASCII): snap2ascii.o
	$(LINK) -o snap2ascii snap2ascii.o
$(SUMM2SNAP): summ2snap.o
	$(LINK) -o summ2snap summ2snap.o

sodapop.o: ./src/evolve.cpp ./src/Gene.h ./src/Cell.h ./src/global.h
	$(COMPILE) -o sodapop.o ./src/evolve.cpp
snap2ascii.o: ./tools/snap2ascii.cpp ./src/global.h
	$(COMPILE) -o snap2ascii.o ./tools/snap2ascii.cpp
summ2snap.o: ./tools/summ2snap.cpp ./src/PolyCell.h 
	$(COMPILE) -o summ2snap.o ./tools/summ2snap.cpp

clean:
	rm -f *.o
