
// mpic++ -std=c++11 -o mpi_heat_equation_solver_with_threshold mpi_heat_equation_solver_with_threshold.cpp 
// mpirun --hostfile hostfile -np 4 mpi_heat_equation_solver_with_threshold


#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "mpi.h"
#define  MASTER   0

using namespace std;

typedef std::chrono::high_resolution_clock Clock;

int heat_equation_cal(int n, int argc, char *argv[]){

  double sum_diff, diff, old_value;
  int numtasks, taskid, task, chunksize, tag1, tag2, tag3, tag4, i, j, row, col, source;
  int flag = 0;
  double threshold =  0.001;

  // Declare input and output array
  double** input_temp = (double **) malloc((n+2) * sizeof(double*));
  double** output_temp = (double **) malloc((n+2) * sizeof(double*));
  // Initiate the 2D array 
  double* in = (double *) malloc((n+2) * (n+2) * sizeof(double));
  double* out = (double *) malloc((n+2) * (n+2) * sizeof(double));
  for(i = 0; i < n+2; ++i){
    input_temp[i] = &(in[(n+2) * i]);  
  }   
  
  for(i = 0; i < n+2; ++i){
    output_temp[i] = &(out[(n+2) * i]);  
  }

  // Initialization
  for (i =0; i < n+2; ++i){
    for (j=0; j < n+2; ++j){
      if(i <= n & j <= n & i > 0 & j > 0){
        if (i ==1 | i == n | j == 1 | j == n){
          input_temp[i][j] = 100;
          output_temp[i][j] = 100;
        }else{
          input_temp[i][j] = -100;
          output_temp[i][j] = -100;
        }
      }else{
        input_temp[i][j] = 0;
        output_temp[i][j] = 0;
      }
    }
  }

  double size = (n+2) * (n+2);
  // MPI stuff
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
  MPI_Status status;

  if (numtasks != 4) {
    printf("Quitting. Number of MPI tasks must be 3.\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
    exit(0);
   }

	// Initiate variable
	double c = 0.1;
	double ds = 1.0 / (n+1);
	double dt = (ds)*(ds) * 1.0 / (4*c);
  chunksize = n * 2 / numtasks;
  tag1 = 1;
  tag2 = 2;
  tag3 = 3;
  tag4 = 4;

  /***** Master task only ******/

  auto t1 = Clock::now();  

  // Iterate till converge
  while(1) { 

    diff = 0;
    sum_diff = 0;
    /* Master does its part of the work */
    if (taskid == MASTER){

      /* Send each task its portion of the array - master keeps 1st part */

      for (int processor = 1; processor <=3; processor++)
      {
        if (processor == 1)
        {
          row = 0;
          col = chunksize;
        }else if (processor == 2)
        {
          row = chunksize;
          col = 0;
        }else if (processor == 3)
        {
          row = chunksize;
          col = chunksize;
        }
        MPI_Send(&row, 1, MPI_INT, processor, tag1, MPI_COMM_WORLD);
        MPI_Send(&col, 1, MPI_INT, processor, tag2, MPI_COMM_WORLD);
        for (i = row; i <= chunksize + 1; i++ )
        {
          MPI_Send(&input_temp[i][col], chunksize + 2, MPI_DOUBLE, processor, tag3, MPI_COMM_WORLD);        
          MPI_Send(&output_temp[i][col], chunksize + 2, MPI_DOUBLE, processor, tag4, MPI_COMM_WORLD);
        } 
      }

      row = 0;
      col = 0;
      // Calculation
      for (i = row + 1; i < row + chunksize+1; ++ i ) {
        for (j = col + 1; j < col + chunksize+1; ++ j) {
          old_value = output_temp[i][j];
          output_temp[i][j] = input_temp[i][j]+c*dt/(ds*ds)*(input_temp[i+1][j]+input_temp[i-1][j]-4*input_temp[i][j]+input_temp[i][j+1]+input_temp[i][j-1]);
          diff += fabs(old_value - output_temp[i][j]);
        }
      }
      swap(input_temp, output_temp);

      /* Wait to receive results from each task */

      for (int processor = 1; processor <=3; processor++)
      {
        MPI_Recv(&row, 1, MPI_INT, processor, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&col, 1, MPI_INT, processor, tag2, MPI_COMM_WORLD, &status);
        for (i = row + 1 ; i < row + chunksize + 1; i++ )
        { 
          MPI_Recv(&input_temp[i][col+1], chunksize, MPI_DOUBLE, processor, tag3, MPI_COMM_WORLD, &status);
          MPI_Recv(&output_temp[i][col+1], chunksize, MPI_DOUBLE, processor, tag4, MPI_COMM_WORLD, &status);    
        } 
      }
      // store the sum of differniation number in sum_diff
      MPI_Reduce(&diff, &sum_diff, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

      if (sum_diff < threshold) {
        printf("Has been converged! sum_diff is %f\n", sum_diff);
        flag = 1;
      }

    } /* end of master */

    /***** Non-master tasks only *****/
    if (taskid > MASTER){

      /* Receive my portion of array from the master task */
      source = MASTER;
      MPI_Recv(&row, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
      MPI_Recv(&col, 1, MPI_INT, source, tag2, MPI_COMM_WORLD, &status);
      for (i = row; i <= chunksize + 1; i++ )
      {
        
        MPI_Recv(&output_temp[i][col], chunksize + 2, MPI_DOUBLE, source, tag4, MPI_COMM_WORLD, &status);
        MPI_Recv(&input_temp[i][col], chunksize + 2, MPI_DOUBLE, source, tag3, MPI_COMM_WORLD, &status);
      } 
      

      // Calculation
      for (i = row + 1; i < row + chunksize + 1; ++ i ) {
        for (j = col + 1; j < col + chunksize + 1; ++ j) {
          old_value = output_temp[i][j];
          output_temp[i][j] = input_temp[i][j]+c*dt/(ds*ds)*(input_temp[i+1][j]+input_temp[i-1][j]-4*input_temp[i][j]+input_temp[i][j+1]+input_temp[i][j-1]);
          diff += fabs(old_value - output_temp[i][j]);
        }
      }
      swap(input_temp, output_temp);

      task = MASTER;
      MPI_Send(&row, 1, MPI_INT, task, tag1, MPI_COMM_WORLD);
      MPI_Send(&col, 1, MPI_INT, task, tag2, MPI_COMM_WORLD);
      for (i = row + 1; i < row + chunksize + 1; i++ )
      {
        
        MPI_Send(&input_temp[i][col+1], chunksize, MPI_DOUBLE, task, tag3, MPI_COMM_WORLD);
        MPI_Send(&output_temp[i][col+1], chunksize, MPI_DOUBLE, task, tag4, MPI_COMM_WORLD);
      } 

      // store the sum of differniation number in sum_diff
      MPI_Reduce(&diff, &sum_diff, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

    } /* end of non-master */

    // broadcast the results to all MPI Ranks
    MPI_Bcast(&flag, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    if (flag)
    {
      break;
    }

  } /* end of iteration */
  auto t2 = Clock::now();


  // Print the time
  if (taskid == MASTER){

      double timer = chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000000.0;
      cout << "N = " << n << "; Processor = " << numtasks << ";  Runtime: "
        << timer
        << " seconds" << std::endl;
  }


  return 0; 
}

int main(int argc, char *argv[]){ 

  MPI_Init(&argc, &argv);
  heat_equation_cal(100, argc, argv);
  MPI_Finalize(); 

	return 0;
}