# Compiler name
F90 = gfortran

# Set your Fortran compiler flags (optimized for gfortran)
F90FLAGS = -fcheck=bounds -Wall -Wextra -O3 

.SUFFIXES:
.SUFFIXES: .exe .o .f90 

.f90.o:
	$(F90) -c $< -o $*.o $(F90FLAGS)  

.o.exe:
	$(F90) $< -o $*.exe $(F90FLAGS)   

rancat_ascii_example.exe: rancat_ascii_example.o rancat_subs.o ran3.o zbrent_subs.o histograms.o rancat_jswml.o user_routines.o cosmology.o
	$(F90) $^ -o rancat_ascii_example $(F90FLAGS)

# Dependencies
rancat_ascii_example.o: rancat_subs.o rancat_jswml.o
rancat_subs.o: cosmology.o user_routines.o histograms.o
rancat_jswml.o: zbrent_subs.o
uniformcat.o: cosmology.o
clusteredcat.o: rancat_subs.o user_routines.o
rancat_.o: rancat_subs.o
rancat2.o: rancat2_subs.o

