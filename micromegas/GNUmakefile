MICROMEGASDIR  = /afs/naf.desy.de/group/atlas/scratch/fittino/theorycodes/NPFitter.mm2.2/micromegas_2.2.CPC.i
#MICROMEGASDIR  = $(HOME)/theorycodes/micromegas_2.4.O

MICROMEGASLIBS22 = $(MICROMEGASDIR)/sources/micromegas.a \
	  	   $(MICROMEGASDIR)/MSSM/lib/aLib.a \
		   $(MICROMEGASDIR)/MSSM/work/work_aux.a \
	   	   $(MICROMEGASDIR)/CalcHEP_src/lib/model_aux.so \
		   $(MICROMEGASDIR)/CalcHEP_src/lib/sqme_aux.so

MICROMEGASLIBS24 = $(MICROMEGASDIR)/sources/micromegas.a \
	  	   $(MICROMEGASDIR)/MSSM/lib/aLib.a \
		   $(MICROMEGASDIR)/MSSM/work/work_aux.a \
                   $(MICROMEGASDIR)/CalcHEP_src/lib/libSLHAplus.a \
                   $(MICROMEGASDIR)/CalcHEP_src/lib/sqme_aux.so \
                   $(MICROMEGASDIR)/CalcHEP_src/lib/dynamic_me.a

MICROMEGASLIBS   = $(MICROMEGASLIBS22)

CXXFLAGS += -I$(MICROMEGASDIR)/MSSM/lib -I$(MICROMEGASDIR)/sources

include $(MICROMEGASDIR)/CalcHEP_src/FlagsForMake
cLib = $(CALCHEP)/lib

all: micromegas.cpp
	$(CXX) $(CXXFLAGS) micromegas.cpp -o micromegas $(MICROMEGASLIBS) $(lDL) -lm
