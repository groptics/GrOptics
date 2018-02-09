#
#VERSION4.0
#30May2016
#
#      C. Duke
#  Grinnell College
#
# dependencies: root (no gsl or tmva)
include Makefile.common

INCLUDEFLAGS  += -I. 



# add INCLUDEFLAGS last 
CXXFLAGS += $(INCLUDEFLAGS)

OBJECTSPTSRC =  VG_Pilot.o
OBJECTS =  VG_Pilot.o GUtilityFuncts.o GDefinition.o

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
		-Y --  *.cpp  *.h \
		-f Makefile.depend0 > /dev/null 2>&1 &&\
	sed "s/^[a-zA-Z0-9]*\///" Makefile.depend0 > Makefile.depend &&\
	rm -f Makefile.depend0

include Makefile.depend

# DO NOT DELETE
#DEPEND LIST DONT DELETE
