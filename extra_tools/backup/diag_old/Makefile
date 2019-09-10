#FC = mpifrtpx
FC = frtpx
FFLAGS = -Kfast,parallel,ocl,optmsg=2 -Qt -Cpp -X9 -fw -fs
FFLAGS_OMP = -Kopenmp

SRC = v17diag.f90
PROG = 'diag.exe'

MYL = ../../lib/
MATH = gkvp_f0.30_math_SSL2

INC = -I/usr/local/fftw/3.3/include
LIB = -SSL2 -L/usr/local/fftw/3.3/lib64 -lfftw3_omp -lfftw3 -lm


diag:	$(SRC) $(MYL)$(MATH).f90
	
	$(FC) $(FFLAGS) $(FFLAGS_OMP) $(SRC) $(MYL)$(MATH).f90 -o $(PROG) $(INC) $(LIB)

	\rm -f *.lst *.o *.mod

clean:
	rm -f *.lst *.o *.mod *.exe go.diag.o* CompileList.txt

clear:
	rm -f ./*.err ./*.out

