include ../make.sys

.SUFFIXES :
.SUFFIXES : .o .f90

solve_rr.x:solver_rr.o
	$(F90) $< -o $@ $(FFLAGS)

.f90.o:
	$(F90) -c $< $(FFLAGS)

clean:
	rm -f *.o *.x *.mod

solve_rr.o:solve_rr.f90

