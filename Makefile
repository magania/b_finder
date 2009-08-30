AA=/home/magania/src/aatrack/AA
DFLAGS=P21
#----------------- aatrack ---------------
#AAMC=/home/magania/Bs/aatrack/P17/AA
#AA17=/home/magania/Bs/aatrack/P17/AA
#AA20=/home/magania/Bs/aatrack/P20/AA
#AA21=/home/magania/Bs/aatrack/P21/AA
#AA=$(AA21)
#DFLAGS=P21

ZLIB_DIR=/lib32

#INCLUDES = -I$(CERNSOURCE_DIR)/cosrc/include -I $(AA)/CLHEP -I $(AA)/src -I $(AA)
AAINCLUDE = -I $(AA)/CLHEP -I $(AA)/src -I $(AA)
AALIBS = -L $(AA)/lib -lAA -lAna -lPhy -lCLHEP -L$(ZLIB_DIR) -lz
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


MYLIBS = obj/DecayMC.o obj/BdJPsiKstarFinder.o obj/BdJPsiKstarMCFinder.o obj/BsJPsiPhiFinder.o obj/BsJPsiPhiMCFinder.o obj/JPsiFinder.o obj/EvtSaver.o obj/PtlSaver.o obj/PhiFinder.o obj/KstarFinder.o obj/VrtSaver.o obj/TagSaver.o obj/PtlFinder.o obj/BhhFinder.o obj/GammaFinder.o obj/UpsilonFinder.o obj/XYGammaFinder.o obj/PiGGFinder.o 
#all: bs_finder bd_finder jpsi_finder hh_finder yp_finder
all: pi_finder;

obj/%.o : src/%.cpp include/%.h
	g++ -m32  $(INCLUDES) -g -o $@ -c $<

$(MYLIBS): obj/%.o : src/%.cpp include/%.h
	g++ -m32 $(INCLUDES) -g -o $@ -c $<

pi_finder : pi_finder.cpp $(MYLIBS) $(AATRACK)
	g++ -m32 $(INCLUDES) -D$(DFLAGS) -g -o obj/pi_finder.o -c $<
	g++ -m32 $(FLAGS) $(INCLUDES) $(LIBS) -g -o pi_finder obj/pi_finder.o $(MYLIBS) $(AATRACK)

yp_finder : yp_finder.cpp $(MYLIBS) $(AATRACK)
	g++ -m32 $(INCLUDES) -D$(DFLAGS) -g -o obj/yp_finder.o -c $<
	g++ -m32 $(FLAGS) $(INCLUDES) $(LIBS) -g -o yp_finder obj/yp_finder.o $(MYLIBS) $(AATRACK)

hh_finder : hh_finder.cpp $(MYLIBS) $(AATRACK)
	g++ -m32 $(INCLUDES) -D$(DFLAGS) -g -o obj/hh_finder.o -c $<
	g++ -m32 $(FLAGS) $(INCLUDES) $(LIBS) -g -o hh_finder obj/hh_finder.o $(MYLIBS) $(AATRACK)
	
bs_finder : bs_finder.cpp $(MYLIBS) $(AATRACK)
	g++ -m32 $(INCLUDES) -D$(DFLAGS) -g -o obj/bs_finder.o -c $<
	g++ -m32 $(FLAGS) $(INCLUDES) $(LIBS) -D$(DFLAGS) -g -o bs_finder_${DFLAGS} obj/bs_finder.o $(MYLIBS) $(AATRACK)

bd_finder : bd_finder.cpp $(MYLIBS) $(AATRACK)
	g++ -m32 $(INCLUDES) -D$(DFLAGS) -g -o obj/bd_finder.o -c $<
	g++ -m32 $(FLAGS) $(INCLUDES) $(LIBS) -D$(DFLAGS) -g -o bd_finder_${DFLAGS} obj/bd_finder.o $(MYLIBS) $(AATRACK)

jpsi_finder : jpsi_finder.cpp $(MYLIBS) $(AATRACK)
	g++ -m32 $(INCLUDES) -D$(DFLAGS) -g -o obj/jpsi_finder.o -c $<
	g++ -m32 $(FLAGS) $(INCLUDES) $(LIBS) -D$(DFLAGS) -g -o jpsi_finder_${DFLAGS} obj/jpsi_finder.o $(MYLIBS) $(AATRACK)

b2mu_ana : b2mu_ana.cpp $(MYLIBS) $(AATRACK)
	g++ -m32 $(INCLUDES) -D$(DFLAGS) -g -o obj/b2mu_ana.o -c $<
	g++ -m32 $(FLAGS) $(INCLUDES) $(LIBS) -D$(DFLAGS) -g -o b2mu_ana obj/b2mu_ana.o $(MYLIBS) $(AATRACK)

#bs_mcfinder : bs_finder.cpp $(MYLIBS) $(AATRACK)
#	g++ $(INCLUDES) -g -DMC -o obj/bs_mcfinder.o -c $<
#	g++ $(FLAGS) $(INCLUDES) $(LIBS) -D$(DFLAGS) -g -o bs_mcfinder obj/bs_mcfinder.o $(MYLIBS) $(AATRACK)


clean: 
	rm $(MYLIBS)

