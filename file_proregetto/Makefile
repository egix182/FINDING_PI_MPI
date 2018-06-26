all: Trapezio MC_Subset MC_Repetition

Trapezio: Trapezoidal_Rule_MPI.c
	mpicc Trapezoidal_Rule_MPI.c -o Trapezio -lm

MC_Subset: Monte_Carlo_Subset_MPI.c
	mpicc Monte_Carlo_Subset_MPI.c -o MC_Subset

MC_Repetition: Monte_Carlo_Repetition_MPI.c
	mpicc Monte_Carlo_Repetition_MPI.c -o MC_Repetition

clean: 
	rm -f  Trapezio MC_Subset MC_Repetition
info: 
	@echo  Run command:
	@echo "mpirun -np [#_of_proc] -host [host_name] [program_name] [some_input_if_required]"
	@echo "mpirun -np [#_of_proc] --hostfile [host_file_name] [program_name] [some_input_if_required]"
