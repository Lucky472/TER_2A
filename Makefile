F90=gfortran
FFLAGS=-O0 -Wall -ffpe-trap=invalid,zero,overflow
EXE=exec

all : $(EXE)

$(EXE) : mod_precision.o mod_caracteristiques.o mod_tri_maillage.o mod_maillage.o mod_resolution.o mod_sortie.o chaleur.o
	$(F90) -o $@ $^

%.o : %.f90
	$(F90) $(FFLAGS) -c $<

clean :
	rm -f *.mod *.o SORTIE/*.vtk $(EXE)