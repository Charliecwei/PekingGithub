/* 
  Foundations of Parallel and Distributed Computing, Falls 2019.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to use MPI send/recv for data exchange. 
*/
#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
  int a, b, i, size, rank;
  float to, t1;

  MPI_Init(&argc, &argv); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  to = MPI_Wtime();
  for (int i = 0; i < 200; ++i)
  {
    if (rank ==0) {
      a = -1;
      MPI_Send(&a, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
      MPI_Recv(&b, 1, MPI_INT, size-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("Process %d received token %d from process %d\n", rank, b, size-1);
    } else if (rank == size-1) {
      MPI_Recv(&b, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("Process %d received token %d from process %d\n", rank, b, rank-1);
      MPI_Send(&b, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); 
    } else {
      MPI_Recv(&b, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("Process %d received token %d from process %d\n", rank, b, rank-1);
      MPI_Send(&b, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }
  }
  t1 = MPI_Wtime();
  t1 = (t1 - to)/(200);
  to = t1/size;
    
if (rank ==0)
{
  printf("The cycle time is %f\n",t1);
  printf("Each tencent time is %f\n ", to);
}

  MPI_Finalize();
  return 0;
}
