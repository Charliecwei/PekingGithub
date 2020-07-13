/* 
  Foundations of Parallel and Distributed Computing, Falls 2019.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to compute matrix multiplicaiton faster. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mat_mul.h"

int main(int argc, char **argv){
  
  int Mwidth = 512;
  int Nwidth = 512;
  int Swidth = 512; // default value
  if (argc > 1) {Mwidth = atoi(argv[1]);
                 Nwidth = atoi(argv[2]);
                 Swidth = atoi(argv[3]);} // user-specified value
  printf("\nMatrix M[%d,%d],N[%d,%d].\n", Mwidth,Nwidth,Nwidth,Swidth);

  srand(time(NULL));
  float* M = rand_mat(Mwidth, Nwidth);
  float* N = rand_mat(Nwidth, Swidth);
  //printf("M %f,\n N %f,%f,%f",M[0],N[0],N[1],N[2]);
  float* cpu_P = raw_mat(Mwidth, Swidth);
  float* gpu_P = raw_mat(Mwidth, Swidth);

  long long cpu_start_time = start_timer();
  cpu_mat_mul(M, N, cpu_P, Mwidth,Nwidth,Swidth);
  long long cpu_time = stop_timer(cpu_start_time, "CPU");

  long long gpu_start_time = start_timer();
  gpu_mat_mul(M, N, gpu_P, Mwidth,Nwidth,Swidth);
  long long gpu_time = stop_timer(gpu_start_time, "GPU");


  // Check the correctness of the GPU results
  int num_wrong = 0;
  for (int i = 0; i < Mwidth * Swidth; i++) {
  //  printf("cup:%f,gpu:%f\n",cpu_P[i],gpu_P[i]);
    if (fabs(cpu_P[i] - gpu_P[i]) > 0.000001) {num_wrong++;
        printf("The wrong number is %d\n",i);}
}	
  // Report the correctness results
  if (num_wrong) printf("GPU %d / %d values incorrect\n", num_wrong, N);
  else           printf("GPU all values correct\n");

}
