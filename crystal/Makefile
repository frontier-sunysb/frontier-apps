#*******************************************************************************
# Valmor de Almeida  dealmeidav@ornl.gov; vfda
# Thu Nov 19 15:07:55 EST 2015; vfda
#*******************************************************************************
#
# There should be no need to change this file
# Don't change unless you know what you are doing
#
XTAL_DIR = .
#
#*******************************************************************************
# Compiler flags and libraries
#*******************************************************************************
#
include options
include local
#
#
#***********************************************************************
#           Form executables
#***********************************************************************
#
#----------------------------------------------------------------------*
# Interpreter component
#----------------------------------------------------------------------*
#
 xtalSrcDir := src
#
 xtalSrcFiles     := $(foreach dir,$(xtalSrcDir),$(wildcard $(dir)/*.cpp))
 xtalSrcFileNames := $(notdir $(xtalSrcFiles))
 xtalObj          := $(addprefix obj/, $(subst .cpp,.o,$(xtalSrcFileNames)))
 xtalHeadFiles    := $(wildcard src/include/*.h)
#
#----------------------------------------------------------------------*
# Targets
#----------------------------------------------------------------------*
#
 crystal : $(xtalObj)
	$(CC) -o crystal.x $(XTAL_NEED_LIB_DIRS) $(CCFLAGS) $(xtalObj) $(XTAL_NEED_LIBS) $(CCLIBS)
#
#----------------------------------------------------------------------*
# Dependencies
#----------------------------------------------------------------------*
#
 $(xtalObj) : $(xtalHeadFiles)
#
#----------------------------------------------------------------------*
#
 VPATH = src/ \
         src/include
#
 obj/%.o : %.cpp
	$(CC) $(XTAL_NEED_INCL_DIRS) $(CCFLAGS) -c $<
	@ mv $(notdir $@) obj/.
#
#***********************************************************************
#           Cleaning targets
#***********************************************************************
#
 clean:
#
	@echo  !!!!!!!!!!!!!!!!!!!
	@echo  !!! CLEANING XTAL!!
	@echo  !!!!!!!!!!!!!!!!!!!
	@/bin/rm -f obj/*.o
	@/bin/rm -f *.x
	@echo  !!!!!!!!!!!!
	@echo  !!! DONE !!!
	@echo  !!!!!!!!!!!!
 debug:
	@echo " "
	@echo "xtalSrcFiles = "
	@echo $(xtalSrcFiles)
	@echo " "
	@echo "xtalSrcFileNames = "
	@echo $(xtalSrcFileNames)
	@echo " "
	@echo "xtalObj = "
	@echo $(xtalObj)
	@echo " "
	@echo "xtalHeadFiles = "
	@echo $(xtalHeadFiles)
