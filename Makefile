# makefile
#VERSION4.0 
#30May2016
#      C. Duke
#  Grinnell College

#
# dependencies: root and gsl
# robast is downloaded, see url in Makefile.common

include Makefile.common

# for ROBAST build (NOTE: ROBAST's include dir is added to INCLUDEFLAGS)
#include Makefile.robast

INCLUDEFLAGS  += -I. -I./include

vpath %.h include
vpath %.cpp src

# add INCLUDEFLAGS  
CXXFLAGS += $(INCLUDEFLAGS)

.PHONY:	all
#.PHONY all
all: robast grReaderFactory grOptics
#all: robast grOptics

# directory to receive all .o files
OBJ := obj

OBJECTS =  $(OBJ)/GUtilityFuncts.o $(OBJ)/GPilot.o \
$(OBJ)/GTelescope.o $(OBJ)/GDCTelescope.o \
$(OBJ)/GTelescopeFactory.o $(OBJ)/GDCTelescopeFactory.o \
$(OBJ)/GReadDCStdBase.o  $(OBJ)/GReadDCStdGrISU.o \
$(OBJ)/GDCGeometry.o  $(OBJ)/GGeometryBase.o \
$(OBJ)/GRayTracerBase.o  \
$(OBJ)/GDCRayTracer.o  \
$(OBJ)/GDefinition.o \
$(OBJ)/GReadPhotonGrISU.o $(OBJ)/GReadPhotonBase.o \
$(OBJ)/GArrayTel.o $(OBJ)/GSimulateOptics.o \
$(OBJ)/GOrderedGrid.o $(OBJ)/GRootDCNavigator.o \
$(OBJ)/GRootWriter.o  $(OBJ)/GSCTelescope.o \
$(OBJ)/GReadSCStd.o  $(OBJ)/GSCTelescopeFactory.o \
$(OBJ)/GSegSCTelescope.o \
$(OBJ)/GSegSCTelescopeFactory.o \
$(OBJ)/GReadSegSCStd.o \
$(OBJ)/GSegmentedMirror.o \
$(OBJ)/GSegmentedObscuration.o

TESTOBJECTS = $(OBJ)/GUtilityFuncts.o $(OBJ)/GDefinition.o 

# set target
testUtilities: $(OBJ)/testUtilities.o  $(TESTOBJECTS)
	@echo "building testUtilities"
	$(LD) $(LDFLAGS) -pg $(LIBS) $^ $(OutPutOpt) $@
	@echo SRCDIR $(SRCDIR)
	@echo INCDIR $(INCDIR)

grOptics: $(OBJ)/grOptics.o $(OBJECTS)
	@echo "building grOptics"
	$(LD) $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@

grReaderFactory: $(OBJ)/grReaderFactory.o $(OBJECTS)
	@echo "building grTestSegReaderFactory"
	$(LD) $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@

%:$(OBJ)/%.o 
	@echo "Building $@ ... "
	$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt) $@
	@echo "Done"
	@echo ""

# rule for any compiling any .cpp file
$(OBJ)/%.o : %.cpp
	@echo "        Compiling $< ... "
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo "Done"

#ifeq ($(ROOTCLING_FOUND),)
# to create root dictionary using rootcling
#else
# to create root dictionary using rootcint
#endif

robast: $(ROOTMAP)

$(ROBAST_TGZ):
	echo Downloading from $(ROBAST_URL) into $(ROBAST_TGZ)
	#curl -L $(ROBAST_URL) > $(ROBAST_TGZ)

$(ROBAST_VER): $(ROBAST_TGZ)
	#tar zxvf $(ROBAST_TGZ)

$(ROOTMAP): 
	cd $(ROBAST_VER);\
	make

var: 
	@echo "ld   $(LD)"
	@echo "archsoflags $(ARCHSOFLAGS)"

includeflags:
	@echo "includeflags: $(INCLUDEFLAGS)"

cxxflags:
	@echo "cxxflags:  $(CXXFLAGS)"

ldflags:
	@echo "ldflags:  $(LDFLAGS)"

libs:
	@echo "libs:  $(LIBS)"

ifeq ($(ROOTCLING),)
  TMPV = 1
endif

testR:
	@echo ROOTCLING_FOUND $(ROOTCLING_FOUND)
	@echo PWD $(PWD)
	@echo TMPV $(TMPV)
	@echo SRCDIR $(SRCDIR)
	@echo INCDIR $(INCDIR)

cleanGrOptics: 
	rm -rf grOptics *.pcm coorTrans src/*Dict* include/*Dict* obj/*.o \
            Makefile.depend 

cleanRobast: 
	rm -rf $(ROBAST_VER) $(ROBAST_TGZ)

clean: cleanGrOptics

DEPEND=echo > Makefile.depend0 &&\
	makedepend -s "\#DEPEND LIST DONT DELETE" -- $(INCLUDEFLAGS) \
		-Y --  src/*.cpp  include/*.h  \
		-f Makefile.depend0 > /dev/null 2>&1 &&\
	sed "s/^[a-zA-Z0-9]*\//obj\//" Makefile.depend0 > Makefile.depend &&\
	rm -f Makefile.depend0

depend:
	$(DEPEND)

Makefile.depend:
	$(DEPEND)

include Makefile.depend

# DO NOT DELETE
#DEPEND LIST DONT DELETE
