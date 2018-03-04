# CIS431-LabAssignment4

2D Heat Equation with OpenMPI CIS 431/531 Parallel Computing

Calculate heat distribution, optimizing by OpenMPI.

For OpenMPI with 4 processors and threshold(mpi_heat_equation_solver_with_threshold.cpp), the compling and running command line is:

~~~
   mpic++ -std=c++11 -o mpi_heat_equation_solver_with_threshold mpi_heat_equation_solver_with_threshold.cpp 
   mpirun --hostfile hostfile -np 4 mpi_heat_equation_solver_with_threshold
~~~

For serial one with 1 processor (heat_equation_solver_1p.cpp), the compling and running command line is:

~~~
   g++ -std=c++11 -o heat_equation_solver_1p heat_equation_solver_1p.cpp
   ./heat_equation_solver_1p
~~~

For OpenMPI with 4 processors and 1000 iterations(mpi_heat_equation_solver_4p.cpp), the compling and running command line is:

~~~
   mpic++ -std=c++11 -o mpi_heat_equation_solver_4p mpi_heat_equation_solver_4p.cpp 
   mpirun --hostfile hostfile -np 4 mpi_heat_equation_solver_4p
~~~

For OpenMPI with 9 processors and 1000 iterations(mpi_heat_equation_solver_9p.cpp), the compling and running command line is:

~~~
   mpic++ -std=c++11 -o mpi_heat_equation_solver_4p mpi_heat_equation_solver_9p.cpp 
   mpirun --hostfile hostfile -np 4 mpi_heat_equation_solver_9p
~~~

I've tested the output, listed below: （unit = second)

~~~
        p = 1   p = 4   p = 9   
        
n=100   0.09    0.17    23.45   

n=1000  10.13   7.48   134.23 

n=5000  226.97  126.43  684.23

~~~

As the data shown above, when n is small(only 100), the cost of data movement outweighs the advantage of parallel(p=4), so serial one won the game. However, as the n gets larger and larger, the speed of openMPI(p=4) is getting faster and faster than serial one. One interesting point is that if p=9, the performance is even much worse than any others. One reason for supporting this evidence is that my Mac doesn't have enough saperate processors to make the parallel work well.
