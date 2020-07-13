#include <stdio.h>
#include <stdlib.h>
void cpu_mat_mul(float* M, float* N, float* P, int Mwidth,int Nwidth,int Swidth) {
  for (int i = 0; i < Mwidth; i++) {
    for (int j = 0; j < Swidth; j++) {
      float sum = 0.0;
      for (int k = 0; k < Nwidth; k++) {
        sum += M[i * Nwidth + k] * N[k * Swidth + j];
      //if (k==2){
        //printf("Row is %d, Col is %d, sum is %.6f\n",i,j, M[i * width + k] * N[k * width + j]);}
      //if (i==2&&j==2){
        //  printf("Row is %d, Col is %d, %.6f,%.6f\n",i,j, M[i * width + k] , N[k * width + j]);}
}


      P[i*Swidth+j] = sum;
      //printf("Row is %d, Col is %d, sum is %.6f\n",i,j, sum);
    }
  }
}
