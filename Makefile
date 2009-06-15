#----------------- aatrack ---------------
AA=/home/magania/src/aatrack/AA
AA17=/home/magania/local/aatrack/aa_p17/AA
AA20=/home/magania/local/aatrack/aa_p20/AA
AA21=/home/magania/local/aatrack/aa_p21/AA
#AA=$(AA21)
DFLAGS=P21

#INCLUDES = -I$(CERNSOURCE_DIR)/cosrc/include -I $(AA)/CLHEP -I $(AA)/src -I $(AA)
AAINCLUDE = -I $(AA)/CLHEP -I $(AA)/src -I $(AA)
AALIBS = -L $(AA)/lib -lAA -lAna -lPhy -lCLHEP -L$(ZLIB_DIR)/lib -lz
AATRACK = $(AA)/lib/libAA.a $(AA)/lib/libAna.a $(AA)/lib/libPhy.a $(AA)/lib/libCLHEP.a

#---------------- root  ------------------
#SPECIALFLAGS=-O -fPIC
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)
ROOTINCLUDE=-I $(shell root-config --incdir)

#RCXX=$(CFLAGS) $(ROOTCFLAGS)
#RLXX=$(LFLAGS) $(ROOTLIBS)

#-----------------------------
INCLUDES = $(AAINCLUDE) $(ROOTINCLUDE) -I include
LIBS = $(AALIBS) $(ROOTLIBS) -lm#-lg2c -lm

#------------------ targets -------------------
FLAGS = $(ROOTCFLAGS) -m32 


MYLIBS = obj/DecayMC.o obj/BsJPsiPhiMCFinder.o obj/JPsiFinder.o obj/EvtSaver.o obj/PtlSaver.o obj/PhiFinder.o obj/VrtSaver.o obj/TagSaver.o 
all: bs_finder

obj/%.o : src/%.cpp include/%.h
	g++ $(INCLUDES) -g -o $@ -c $<

$(MYLIBS): obj/%.o : src/%.cpp include/%.h
	g++ $(INCLUDES) -g -o $@ -c $<
	
bs_finder : bs_finder.cpp $(MYLIBS) $(AATRACK)
	g++ $(INCLUDES) -g -o obj/bs_finder.o -c $<
	g++ $(FLAGS) $(INCLUDES) $(LIBS) -D$(DFLAGS) -g -o bs_finder obj/bs_finder.o $(MYLIBS) $(AATRACK)

clean: 
	rm $(MYLIBS) obj/b_tag.o b_tag

