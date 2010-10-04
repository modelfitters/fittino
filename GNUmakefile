###########################################################################
#                                GNUmakefile
#                            -------------------    
#   Makefile for the Fittino program
#                            -------------------
#   $Id$
#   $Name$
#                            -------------------
#   begin                : September 14, 2003
#   authors              : Philip Bechtle, Peter Wienemann
#   email                : philip.bechtle@desy.de, peter.wienemann@desy.de
###########################################################################

###########################################################################
#                                                                         #
#   This program is free software; you can redistribute it and/or modify  #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation; either version 2 of the License, or     #
#   (at your option) any later version.                                   #
#                                                                         #
###########################################################################

ARCH              = $(shell uname)

HBDIR             = $(HOME)/theorycodes/HiggsBounds-1.2.0/HiggsBounds-f77
#HBDIR             = /afs/naf.desy.de/user/p/prudent/scratch/HiggsBounds-1.2.0/HiggsBounds-f77

ROOTCFLAGS        = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS          = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTLIBS         += -lMinuit
ROOTGLIBS         = $(shell $(ROOTSYS)/bin/root-config --glibs)
LIBS              = $(ROOTLIBS)


CXX               = g++
CXXFLAGS          = -g -Wall -pedantic -Wno-long-long \
                    -Wshadow -fPIC -I. $(ROOTCFLAGS) 

HBLIB		  = -L$(HBDIR) -lHB
LIBGFORTRAN       = -lgfortran
#LIBS	    	 += $(HBLIB) $(LIBGFORTRAN)
#CXXFLAGS         += -DUSELIBHB

#HBLIB             = -L/afs/naf.desy.de/group/atlas/scratch/fittino/theorycodes/HiggsBounds-1.2.0/HiggsBounds-f77/ -lHB
#LIBGFORTRAN       = -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/ -lgfortran
#LIBGFORTRAN      = /usr/lib/libgfortran.so.1
LIBS             += $(HBLIB) $(LIBGFORTRAN)
CXXFLAGS         += -DUSELIBHB

LD                = g++
LDFLAGS           = -g -Wall -pedantic -Wno-long-long -Wshadow
# -pedantic -ansi

SOURCES           = $(wildcard *.cpp)
OBJECTS_BASENAME  = $(SOURCES:.cpp=.o)
OBJECTS_BASENAME += lex.yy.o y.tab.o
INCLUDES          = $(wildcard *.h)
DEPENDS_BASENAME  = $(SOURCES:.cpp=.d)

TARGET            = fittino

all: 	$(TARGET)

## Include dependency files
-include $(SOURCES:.cpp=.d)

%.d: %.cpp
	@echo Generating dependency file ---------- $@ ----------
	@$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) $< \
	              | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
	              [ -s $@ ] || rm -f $@'
#	-rm -f $@
#	$(CXX) -MM $(CXXFLAGS) $< > $@ 

fittino: $(OBJECTS_BASENAME) EndpointFormulas.o
	$(CXX) -o fittino $(LDFLAGS) $(OBJECTS_BASENAME) EndpointFormulas.o $(LIBS)


%.o: %.cpp 
	$(CXX) -c -o $@ $< $(CXXFLAGS)

%.o: %.C
	$(CXX) -c -o $@ $< $(CXXFLAGS)

lex.yy.o: y.tab.c lex.yy.c
	$(CXX) -I. -c -o lex.yy.o lex.yy.c

y.tab.o:  y.tab.c lex.yy.c
	$(CXX) $(CXXFLAGS) -c -o y.tab.o y.tab.c

lex.yy.c: readvalues.lex y.tab.c
	flex -i readvalues.lex 
y.tab.c: readvalues.yacc 
	bison -y -d readvalues.yacc

libsusygen_call_test.a: susygen_call_test.f
	fort susygen_call_test.f
	ar -ruc libsusygen_call_test.a susygen_call_test.o

clean:
	-rm -f *.o *.d *~ lex.yy.c y.tab.h y.tab.c fittino
