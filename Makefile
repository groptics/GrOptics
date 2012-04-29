# makefile
#VERSION2.1 
#1March2012
#      C. Duke
#  Grinnell College
#   Version 1.0
#   June 22, 2010
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
all: robast grOptics

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
           $(OBJ)/GRootDCNavigatorDict.o \
           $(OBJ)/GRootWriterDict.o

# set target

grOptics: $(OBJ)/grOptics.o $(OBJECTS)
	@echo "building grOptics"
	$(LD) $(LDFLAGS) -pg $(LIBS) $^ $(OutPutOpt) $@

%:$(OBJ)/%.o 
	@echo "Building $@ ... "
	$(LD) $(LDFLAGS) -pg  $^ $(LIBS) $(OutPutOpt) $@
	@echo "Done"
	@echo ""

# rule for any compiling any .cpp file
$(OBJ)/%.o : %.cpp
	@echo "        Compiling $< ... "
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo "Done"

# to create root dictionary using rootcint
src/GRootDCNavigatorDict.cpp: GRootDCNavigator.h GRootDCNavigatorLinkDef.h
	@echo "Generating dictionary $< ... "
	rootcint -v -f $@ -c -p $^
	mv src/GRootDCNavigatorDict.h include/.
	@echo "Done"

src/GRootWriterDict.cpp: GRootWriter.h GRootWriterLinkDef.h
	@echo "Generating dictionary $< ... "
	rootcint -v -f $@ -c -p $^
	mv src/GRootWriterDict.h include/.
	@echo "Done"



robast: $(ROOTMAP)

$(ROBAST_TGZ):
	echo Downloading $(TGZ)
	curl -L $(ROBAST_URL) > $(ROBAST_TGZ)

$(ROBAST_VER): $(ROBAST_TGZ)
	tar zxvf $(ROBAST_TGZ)

$(ROOTMAP): $(ROBAST_VER)
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

cleanGrOptics: 
	rm -rf grOptics coorTrans src/*Dict*.cpp include/*Dict*.h obj/*.o \
            Makefile.depend 

cleanRobast: 
	rm -rf $(ROBAST_VER) $(ROBAST_TGZ)

clean: cleanGrOptics cleanRobast

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
