# makefile
#VERSION4.0 
#30May2016
#      C. Duke
#  Grinnell College

#
# dependencies: root and gsl
# robast is downloaded, see url in Makefile.common

include Makefile.common

INCLUDEFLAGS  += -I. -I./include

# directory to receive all .o files
OBJDIR := obj

SRCDIR := src
SrcSuf := cpp
ObjSuf := o

vpath %.h include
vpath %.cpp src

# add INCLUDEFLAGS  
CXXFLAGS += $(INCLUDEFLAGS)

.PHONY:	all
all: robast grReaderFactory grOptics

# make list of all G* source files
SRCS    :=      $(filter-out $(SRCDIR)/GArray_Tel_noGraphs_4.0deg.%,$(wildcard $(SRCDIR)/G*.$(SrcSuf)))
# make list of object files to create, DICTO is currently empty
OBJECTS  := $(patsubst $(SRCDIR)/%.$(SrcSuf),$(OBJDIR)/%.$(ObjSuf),$(SRCS))  $(DICTO) 

#------------------------------------------------------------------
.PHONY: print
print:
	@echo SRCS $(SRCS)
	@echo OBJECTS $(OBJECTS)
	@echo OBJDIR $(OBJDIR)

TESTOBJECTS = $(OBJDIR)/GUtilityFuncts.o $(OBJDIR)/GDefinition.o 

# set target
testUtilities: $(OBJDIR)/testUtilities.o  $(TESTOBJECTS)
	@echo "building testUtilities"
	$(LD) $(LDFLAGS) -pg $(LIBS) $^ $(OutPutOpt) $@
	@echo SRCDIR $(SRCDIR)
	@echo INCDIR $(INCDIR)

grOptics: $(OBJDIR)/grOptics.o $(OBJECTS)
	@echo "building grOptics"
	$(LD) $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@

grReaderFactory: $(OBJDIR)/grReaderFactory.o $(OBJECTS)
	@echo "building grTestSegReaderFactory"
	$(LD) $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@

%:$(OBJDIR)/%.o 
	@echo "Building $@ ... "
	$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt) $@
	@echo "Done"
	@echo ""

# rule for any compiling any .cpp file
$(OBJDIR)/%.o : %.cpp
	@echo "        Compiling $< ... "
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo "Done"

#ifeq ($(ROOTCLING_FOUND),)
# to create root dictionary using rootcling
#else
# to create root dictionary using rootcint
#endif

.PHONY: robast
robast:
ifeq ($(filter $(ROBAST_DIR), $(wildcard */)),)
	@echo download ROBAST tar file, no $(ROBAST_VER) directory
	curl -L $(ROBAST_URL) > $(ROBAST_TGZ)
	tar zxvf $(ROBAST_TGZ)
	@echo copying ROBASTEditedFile/AOpticalComponent.cxx to robast src directory
	cp ROBASTEditedFile/AOpticalComponent.cxx ROBAST-2.4.1/src/.
else
	@echo  $(ROBAST_DIR) directory already exists
endif
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
	rm -rf grOptics *.pcm coorTrans $(SRCDIR)/*Dict* $(INCDIR)/*Dict* $(OBJDIR)/*.o \
            Makefile.depend 

cleanRobast: 
	cd $(ROBAST_VER); make clean

removeRobast:
	rm -rf $(ROBAST_VER) $(ROBAST_TGZ)

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
