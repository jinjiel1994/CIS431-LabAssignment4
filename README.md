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
        
n=100   0.09    1.05    1.06     

n=1000  10.13  81.44   88.38  

n=5000  226.97   6853.9  6922.5 

~~~

I cannot run 1e5*1e5 in my Mac because the magnitude is too large to compute. However, from the gragh listed above, we can see when the magnitude is not quite large, optimization is not significant, while if the magnitude is large enough, multiple threads can run much faster.
