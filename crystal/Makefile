CXX=mpicxx 
F77=mpif77 
F77_LIBS=

libext = 
incs =   -I../../frontier -I../../frontier/src/util   -I/usr/local/ompi/include -D__MPI__   -DHAS_FENV    -D__GD__ -I/usr/lib64/include  -D__HYPRE__  
libincs =  -L/usr/lib64/lib  -L/usr/local/ompi/lib  -L/usr/local/petsc/lib -L/usr/local/hdf4/lib
libs = -lgd   -lmpi -lpthread -lmfhdf -ldf  -ljpeg -lz
build_cpu = x86_64
incsloc = -I. -I../solver -I../../frontier/include

PETSC = yes
PETSC_ARCH=
PETSC_DIR=

Petsc_Include= -I/usr/local/petsc/include
Petsc_Lib= -lHYPRE -lpetsc
crystal : 

solver.o : ../solver/solver.h ../solver/solver.cpp
	${CXX} -c ../solver/solver.cpp $(incsloc) $(incs) ${Petsc_Include}

parab.o : ../solver/solver.h ../solver/parab.cpp
	${CXX} -c ../solver/parab.cpp $(incsloc) $(incs) ${Petsc_Include}

ccartsn.o : crystal.h ../solver/solver.h ccartsn.cpp
	${CXX} -c ccartsn.cpp $(incsloc) $(incs) ${Petsc_Include}

crysub.o : crystal.h ../solver/solver.h crysub.cpp
	${CXX} -c crysub.cpp $(incsloc) $(incs) ${Petsc_Include}

crystal: crystal.h crystal_state.h crystal.cpp solver.o parab.o ccartsn.o crysub.o
	${CXX} -c crystal.cpp -I../../frontier/include $(incs) ${Petsc_Include}
	$(CXX) -o crystal crystal.o -I../../frontier/include -L../../frontier/lib/${build_cpu} solver.o parab.o ccartsn.o crysub.o -L../../frontier/lib/${build_cpu} -lFronTier $(libincs) $(libs) ${Petsc_Lib} ${F77_LIBS}

-include ../../frontier/devel-deps.inc

clean:
	rm -rf *.o crystal 
output :
	rm -rf out/*
distclean:
	rm -rf *.o crystal out/*
cleanall:
	rm -rf *.o crystal out/*
gif :
	hdf2gif out/intfc.hdf out/intfc.gif
TagsFile:
	ctags *.h *.cpp ../solver/solver.h ../solver/solver.cpp ../solver/parab.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]
tagsfile:
	ctags *.h *.cpp ../solver/solver.h ../solver/solver.cpp ../solver/parab.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]
etagsfile:
	etags *.h *.cpp ../solver/solver.h ../solver/solver.cpp ../solver/parab.cpp ../src/front/*.[ch] ../src/intfc/*.[ch] ../src/util/*.[ch]
