#*******************************************************************************
# Valmor de Almeida  dealmeidav@ornl.gov; vfda
# Thu Nov 19 15:07:55 EST 2015; vfda
#*******************************************************************************
#
# There should be no need to change this file
# Don't change unless you know what you are doing
#
SOLV_DIR = .
#
#*******************************************************************************
# Compiler flags and libraries
#*******************************************************************************
#
include options
include local
#
#***********************************************************************
#           Form executables
#***********************************************************************
#
#----------------------------------------------------------------------*
# Interpreter component
#----------------------------------------------------------------------*
#
 solvSrcDir := src
#
 solvSrcFiles     := $(foreach dir,$(solvSrcDir),$(wildcard $(dir)/*.cpp))
 solvSrcFileNames := $(notdir $(solvSrcFiles))
 solvObj          := $(addprefix obj/, $(subst .cpp,.o,$(solvSrcFileNames)))
 solvHeadFiles    := $(wildcard src/include/*.h)
#
#----------------------------------------------------------------------*
# Targets
#----------------------------------------------------------------------*
#
 solvlib : $(solvObj)
	ar rv libsolver.a $(SOLV_NEED_LIB_DIR) $(SOLV_NEED_LIB) $(solvObj)
	ranlib libsolver.a
	mv libsolver.a lib/.
#
#----------------------------------------------------------------------*
# Dependencies
#----------------------------------------------------------------------*
#
 $(solvObj) : $(solvHeadFiles)
#
#----------------------------------------------------------------------*
#
 VPATH = src/ \
         src/include
#
 obj/%.o : %.cpp
	$(CC) $(SOLV_NEED_INCL_DIRS) $(CCFLAGS) -c $<
	@mv $(notdir $@) obj/.
#
#***********************************************************************
#           Cleaning targets
#***********************************************************************
#
 clean:
#
	@echo  !!!!!!!!!!!!!!!!!!!
	@echo  !!! CLEANING SOLV!!
	@echo  !!!!!!!!!!!!!!!!!!!
	@/bin/rm -f obj/*.o
	@/bin/rm -f lib/lib*.a
	@echo  !!!!!!!!!!!!
	@echo  !!! DONE !!!
	@echo  !!!!!!!!!!!!
 debug:
	@echo " "
	@echo "solvSrcFiles = "
	@echo $(solvSrcFiles)
	@echo " "
	@echo "solvSrcFileNames = "
	@echo $(solvSrcFileNames)
	@echo " "
	@echo "solvObj = "
	@echo $(solvObj)
	@echo " "
	@echo "solvHeadFiles = "
	@echo $(solvHeadFiles)

