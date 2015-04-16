#  version: 6.1.0    date: 2014-05-06
# template makefile
#      C. Duke
#  Grinnell College
#   Version 2.0
#   May, 2011
#
# dependencies: root (no gsl or tmva)
include Makefile.common

INCLUDEFLAGS  += -I. 



# add INCLUDEFLAGS last 
CXXFLAGS += $(INCLUDEFLAGS)

OBJECTSPTSRC =  VG_Pilot.o
OBJECTS =  VG_Pilot.o GUtilityFuncts.o


# define directory to receive all .o files
#OBJ := obj

# pattern rule for any application code, e.g. vg_get_data
# where a corresponding .cpp code must exist in a VPATH directory
# The rule updates the grisep library if necessary
# links the application code to library routines
# and removes the .o application file.

makePhotonsOffset:: makePhotonsOffset.o  $(OBJECTS)
	@echo "linking makePhotonsOffset"
	$(LD) $(LDFLAGS) -g  $^  $(LIBS)  $(OutPutOpt) $@

makePhotonsPtSrc: makePhotonsPtSrc.o $(OBJECTSPTSRC)
	@echo "linking makePhotonsPtSrc"
	$(LD) $(LDFLAGS) -g  $^  $(LIBS)  $(OutPutOpt) $@

%:%.o 
	@printf "Building $@ ... "
	$(LD) $(LDFLAGS) -g  $^ $(LIBS) $(OutPutOpt) $@
	@echo "Done"
	@echo ""


# rule for any compiling any .cpp file
%.o : %.cpp  
	@printf "Compiling $< ... "
	@$(CXX) $(CXXFLAGS) -c $<
	@echo "Done"

#VGA_GeoTelescopeDict.cpp: VGA_GeoTelescope.h VGA_GeoTelescopeLinkDef.h
#	@printf "Generating dictionary $< ... "
#	rootcint -v -f $@ -c -p $^
#	@echo "Done"

var: 
	@echo "ld   $(LD)"
	@echo "archsoflags $(ARCHSOFLAGS)"

cxxflags:
	@echo "cxxflags:  $(CXXFLAGS)"

ldflags:
	@echo "ldflags:  $(LDFLAGS)"

libs:
	@echo "libs:  $(LIBS)"


clean:
	rm -f *.o makePhotonsPtSrc makePhotonsOffset Makefile.depend 
	rm -f Makefile.depend0.bak

depend:
	$(DEPEND)

Makefile.depend:
	$(DEPEND)

DEPEND=echo > Makefile.depend0 &&\
	makedepend -s "\#DEPEND LIST DONT DELETE" -- $(INCLUDEFLAGS) \
		-Y --  src/*.cpp  include/*.h \
		-f Makefile.depend0 > /dev/null 2>&1 &&\
	sed "s/^[a-zA-Z0-9]*\///" Makefile.depend0 > Makefile.depend &&\
	rm -f Makefile.depend0

include Makefile.depend

# DO NOT DELETE
#DEPEND LIST DONT DELETE
