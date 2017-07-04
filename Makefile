RM =rm
CP =cp
CC =g++

IDIR =./include
CXXFLAGS =-std=c++11 -Wall -O2 -I$(IDIR)
EXEC =sodapop
LINK = $(CXX) $(CXXFLAGS)
COMPILE = $(CXX) $(LIBS) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -c

SODAPOP = soda_pop
SNAP2ASCII = soda_snap
SUMM2SNAP = soda_summ

INSTALLDIR = /usr/local/bin

all: $(SODAPOP) $(SNAP2ASCII) $(SUMM2SNAP)
install:
	$(CP) $(SODAPOP) $(INSTALLDIR)/
	$(CP) $(SNAP2ASCII) $(INSTALLDIR)/
	$(CP) $(SUMM2SNAP) $(INSTALLDIR)/
	
uninstall:
	$(RM) -r $(INSTALLDIR)/$(SODAPOP)
	$(RM) -r $(INSTALLDIR)/$(SNAP2ASCII)
	$(RM) -r $(INSTALLDIR)/$(SUMM2SNAP)

$(SODAPOP): sodapop.o
	$(LINK) -o soda_pop sodapop.o
$(SNAP2ASCII): snap2ascii.o
	$(LINK) -o soda_snap snap2ascii.o
$(SUMM2SNAP): summ2snap.o
	$(LINK) -o soda_summ summ2snap.o

sodapop.o: ./src/evolve.cpp ./src/Gene.h ./src/Cell.h ./src/global.h
	$(COMPILE) -o sodapop.o ./src/evolve.cpp
snap2ascii.o: ./tools/snap2ascii.cpp ./src/global.h
	$(COMPILE) -o snap2ascii.o ./tools/snap2ascii.cpp
summ2snap.o: ./tools/summ2snap.cpp ./src/PolyCell.h 
	$(COMPILE) -o summ2snap.o ./tools/summ2snap.cpp

clean:
	rm -f *.o
