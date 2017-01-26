# makefile
#VERSION4.0 
#30May2016
#      C. Duke
#  Grinnell College
# UPDATED VERSION
# Makefile for grOptics.  requires root and robast. robast downloaded
# using curl if necessary.  Code changes to robast classes necessary for
# using robast with GrOptics copied to robast prior to robast compilation
# (see directory GrOptics/ROBASTEditedFile contents)
#
# dependencies: root and gsl
# robast is downloaded, see url in Makefile.common

include Makefile.common

INCLUDEFLAGS  += -I. -I./include
# add INCLUDEFLAGS   
CXXFLAGS += $(INCLUDEFLAGS)

# directory to receive all .o files
OBJDIR := obj
SRCDIR := src
INCDIR := include
SrcSuf := cpp
ObjSuf := o

DICTO = $(OBJDIR)/GrOpticsDict.o 
DICTS = $(SRCDIR)/GrOpticsDict.cpp 
DICTI = $(INCDIR)/GrOpticsDict.h 

# make listing of all source files, add dictionary file will generate later
FILTEROUT =  $(SRCDIR)/GArray_Tel_noGraphs_4.0deg.cpp $(SRCDIR)/GrOpticsDict.cpp  
SRCS    :=      $(filter-out $(FILTEROUT),$(wildcard $(SRCDIR)/G*.$(SrcSuf))) $(DICTS)

# make list of object files to create, DICTO is currently empty
OBJS  := $(patsubst $(SRCDIR)/%.$(SrcSuf),$(OBJDIR)/%.$(ObjSuf),$(SRCS))

#make listing of include files
INCS := $(filter-out $(INCDIR)LinkDef.h, $(wildcard $(INCDIR)/G*.h ) ) 

# make listing of include files for dictionary creation
INCD = $(INCDIR)/GRootWriter.h 

.PHONY: all printDebug robast var includeflags cxxflags ldflags \
	libs cleanGrOptics cleanRobast removeRobast testDict

all: robast grOptics

estDict: $(DICTS)
	@echo in testDict

#------------------------------------------------------------------
printDebug:
	@echo SRCS $(SRCS)
	@echo OBJS $(OBJS)
	@echo OutPutOpt $(OutPutOpt)
	@echo ROOTCLING_FOUND $(ROOTCLING_FOUND)
	@echo PWD $(PWD)
	@echo SRCDIR $(SRCDIR)
	@echo INCDIR $(INCDIR)
	@echo OBJDIR $(OBJDIR)

TESTOBJECTS = $(OBJDIR)/GUtilityFuncts.o $(OBJDIR)/GDefinition.o 

testUtilities: $(OBJDIR)/testUtilities.o  $(TESTOBJECTS)
	@echo "building testUtilities"
	$(LD) $(LDFLAGS) -pg $(LIBS) $^ $(OutPutOpt) $@
	@echo SRCDIR $(SRCDIR)
	@echo INCDIR $(INCDIR)

grOptics: $(OBJDIR)/grOptics.o $(OBJS)
	@echo "building grOptics"
	$(LD) $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@

grReaderFactory: $(OBJDIR)/grReaderFactory.o $(OBJS)
	@echo "building grTestSegReaderFactory"
	$(LD) $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@


# rule for any compiling any .cpp file
$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@echo "        Compiling $< ... "
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo "Done"


ifneq ($(ROOTCLING_FOUND),)
# to create root dictionary using rootcling
$(DICTS): $(INCD) $(INCDIR)/LinkDef.h
	@echo "Generating dictionary with rootcling $< ... "
	rootcling -v -f $@ -c  $^
	@cp src/GrOpticsDict_rdict.pcm .	
	@echo "Done"

else
# to create root dictionary using rootcint
$(DICTS): $(INCD) $(INCDIR)/LinkDef.h
	@echo "Generating GrOpticsDict dictionary $< ... "
	rootcint -v -f $@ -c  $^
	@echo "rootcint SHELL mv src/GRootWriterDict.h include/."
	cp src/GrOpticsDict.h  include/.
	@echo "Done"
endif

# -------------  ROBAST -----------------------------------
# search for ROBAST_VER directory

# get edit files for ROBAST from robast edit directory
EDIT_DIR := ROBASTEditedFile
INCFILE := $(wildcard $(EDIT_DIR)/*.h)
CXXFILE := $(wildcard $(EDIT_DIR)/*.cxx)
COPY_FILE = $(EDIT_DIR)/copyFile

# see if ROBAST is already installed, look for directory
ROBDIR = $(wildcard $(ROBAST_VER))
#ROBDIR = $(strip $(wildcard $(ROBAST_VER))) 

ifeq ($(ROBDIR),)
robast:
	@echo did not find $(ROBAST_VER) directory
	@echo     ready to download ROBAST into tar file
	curl -L $(ROBAST_URL) > $(ROBAST_TGZ)
	@echo     ready to untar ROBAST into $(ROBAST_VER) directory
	tar zxf $(ROBAST_TGZ)
	@echo     ready to copy robast edit files
	$(foreach f,$(CXXFILE), cp $(f) $(ROBAST_VER)/src/.; )
	$(foreach f,$(INCFILE), cp $(f) $(ROBAST_VER)/include/.;) 
	touch $(COPY_FILE)
	cd $(ROBAST_VER); make
	@echo    robast ready and waiting

else
.PHONY: robast
robast: $(COPY_FILE)
	@echo    robast ready and waiting
endif

$(COPY_FILE): $(CXXFILE) $(INCFILE)
	@echo     ready to copy robast edit files
	$(foreach f,$(filter %.cxx, $?), cp $(f) $(ROBAST_VER)/src/.;) 
	$(foreach f,$(filter %.h, $?), cp $(f) $(ROBAST_VER)/include/.;) 

	touch $(COPY_FILE)
	cd $(ROBAST_VER); make;

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

objs:
	@echo $(OBJS)

srcs:
	@echo $(SRCS)
incs:
	@echo $(INCS)

incd:
	@echo $(INCD)
dicts:
	@echo $(DICTS)
dicti:
	@echo $(DICTI)

cleanGrOptics: 
	rm -rf grOptics *.pcm $(SRCDIR)/*Dict* $(INCDIR)/*Dict* $(OBJDIR)/*.o \
            Makefile.depend grReaderFactory $(COPY_FILE) *.pcm

cleanRobast:
	cd $(ROBAST_VER); make clean;
	rm -f  $(COPY_FILE)

removeRobast:
	rm -rf $(ROBAST_VER) $(ROBAST_TGZ)
	rm -f  $(COPY_FILE)

clean: cleanGrOptics

DEPEND=echo > Makefile.depend0 &&\
	makedepend -s "\#DEPEND LIST DONT DELETE" -- $(INCLUDEFLAGS) \
		-Y --  $(SRCDIR)/*.cpp  $(INCDIR)/*.h  \
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


