F90=gfortran
FFLAGS=-O2
EXE=exec

all : $(EXE)

$(EXE) : mod_precision.o mod_tri_maillage.o mod_maillage.o mod_sortie.o chaleur.o
	$(F90) -o $@ $^

%.o : %.f90
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f *.mod *.o SORTIE/*.vtk $(EXE)