RM =rm
CP =cp
CC =g++

IDIR =./include
CXXFLAGS =-std=c++11 -Wall -O3 -I$(IDIR)
LINK = $(CXX) $(CXXFLAGS)
COMPILE = $(CXX) $(LIBS) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -c

SODAPOP = sodapop
SNAP2ASCII = sodasnap
SUMM2SNAP = sodasumm

LZ4_DIR = src/lz4-dev/
INSTALLDIR = /usr/local/bin

.PHONY: lz4_code $(LZ4_DIR)

all: $(SODAPOP) $(SNAP2ASCII) $(SUMM2SNAP) lz4_code

install:
	@echo \#\#\# Installing binaries to $(INSTALLDIR)/...
	$(CP) $(SODAPOP) $(INSTALLDIR)/
	$(CP) $(SNAP2ASCII) $(INSTALLDIR)/
	$(CP) $(SUMM2SNAP) $(INSTALLDIR)/
	
uninstall:
	@echo \#\#\# Uninstalling binaries from $(INSTALLDIR)/...
	$(RM) -r $(INSTALLDIR)/$(SODAPOP)
	$(RM) -r $(INSTALLDIR)/$(SNAP2ASCII)
	$(RM) -r $(INSTALLDIR)/$(SUMM2SNAP)

$(SODAPOP): sodapop.o rng.o
	$(LINK) -o sodapop sodapop.o rng.o
$(SNAP2ASCII): snap2ascii.o
	$(LINK) -o sodasnap snap2ascii.o rng.o
$(SUMM2SNAP): summ2snap.o
	$(LINK) -o sodasumm summ2snap.o rng.o
$(LZ4_DIR): 
	make -C $@

lz4_code: $(LZ4_DIR)

rng.o: ./src/rng.cpp
	$(COMPILE) -o rng.o ./src/rng.cpp
sodapop.o: ./src/evolve.cpp ./src/Gene.h ./src/Cell.h ./src/global.h
	$(COMPILE) -o sodapop.o ./src/evolve.cpp
snap2ascii.o: ./tools/snap2ascii.cpp ./src/global.h
	$(COMPILE) -o snap2ascii.o ./tools/snap2ascii.cpp
summ2snap.o: ./tools/summ2snap.cpp ./src/PolyCell.h 
	$(COMPILE) -o summ2snap.o ./tools/summ2snap.cpp

clean:
	rm -f *.o
	make -C $(LZ4_DIR) clean
