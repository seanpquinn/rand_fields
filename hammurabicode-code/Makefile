#
#  This is the hammurabi Makefile.  Since I use multiple systems with
#  multiple makefiles, I split them.  In these files are the things
#  that change.  The rest shouldn't need to be any different.

include Makefile.custom.linux
#include Makefile.custom.servso1
#include Makefile.custom.magique3



#------------------------------------------------------------
#------------------------------------------------------------
#From here onwards, nothing should need to be changed 
#------------------------------------------------------------
#------------------------------------------------------------

CXXFLAGS+= $(ALL_INC) -c -DGALDEF_PATH=\"./GALDEF\" -DFITSDATA_PATH=\"./FITS/\" -DDATA_PATH=\"./DATA\" 
LFLAGS = $(ALL_L) -lhammurabi -lNE2001 -lhealpix_cxx -lcxxsupport -lc_utils -lfftpack -lcfitsio -lgsl -lgslcblas -lm -lfftw3 -lgfortran $(LFLAGS_CUSTOM)

all: hammurabi print_constants


%.o : %.cpp
	$(CXX) $(CXXFLAGS)  -o $@ $<
%.o : %.cc
	$(CXX) $(CXXFLAGS) -o $@ $<
%.o : %.c
	$(CXX) $(CXXFLAGS) -o $@ $<	
%.o : %.f
	$(FC) $(FCFLAGS) -o $@ $<

class_TE_density.o: class_TE_density.cpp
	$(CXX) $(CXXFLAGS) $(LFORTFLAGS) -o $@ $<

NE2001_OBJ = dmdsm.NE2001.o density.NE2001.o neclumpN.o nevoidN.o neLISM.NE2001.o

HAMMURABI_OBJ = $(NE2001_OBJ) class_List.o class_Integrator.o class_B_field2.o tess_tools.o class_TE_density.o class_CRE.o namespace_Vec_Handling.o

run:
	mkdir run

$(HAMMURABI_OBJ):	hammurabi.h CGS_units_file.h proto*.h tess_tools.h

libNE2001.a: $(NE2001_OBJ)
	ar rc $@ $(NE2001_OBJ)
	ranlib $@

hammurabi: 	$(HAMMURABI_OBJ) $(NE2001_DATA) libNE2001.a libhammurabi.a run hammurabi.o
	$(CXX) -o run/hammurabi hammurabi.o $(LFLAGS) $(LFORTFLAGS)

print_constants: print_constants.o
		 $(CXX) -o run/print_constants print_constants.o 

# To link in B_field to galprop, use this:
libhamprop.a:  class_B_field2.o tess_tools.o namespace_Vec_Handling.o 
	ar cru libhamprop.a class_B_field2.o tess_tools.o namespace_Vec_Handling.o 

libhammurabi.a: $(HAMMURABI_OBJ)
	ar cru libhammurabi.a $(HAMMURABI_OBJ)


clean:
	rm *.a *.o run/hammurabi

test:
	(cd unit_test; ./test.csh)


tarfile:
	tar cvzf hammurabi_mine.tgz *cpp *h Makefile Makefile.custom.mac Makefile.custom.servso1 *.f 
