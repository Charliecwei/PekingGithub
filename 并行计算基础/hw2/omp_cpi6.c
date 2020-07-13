/* 
 Foundations of Parallel and Distributed Computing, Falls 2019.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to calculate pi with parallel construct. 
*/

#include <omp.h>
#include <stdio.h>
#include <math.h>

#define PI25DT 3.141592653589793238462643
#define TOL 1e-16
#define SS 0.216506350946110
struct Sd{int value;double value2;};
struct Sd adapquad(double a, double b);
double trap(double a, double b);
int main(int argc, char *argv[]){
  int    nthreads, n, i;
  double pi,  a, b, t0, t1;
  struct Sd result;
  t0 = omp_get_wtime();
  pi = 0.0;
  n = 0;
#pragma omp parallel default(shared) \
  private(i, a, b, result) reduction(+:pi) \
  reduction(+:n)
  {
    nthreads = omp_get_num_threads();
    double h = 0.5/(double)nthreads; 
#pragma omp for schedule(static,1) 
    for (i = 1; i <= nthreads; i++) {
      a = h*((double)(i-1));
      b = h*((double)i);
      result = adapquad(a,b);
      pi = result.value2;
      n = result.value;
      
    }
  }
  t1 = omp_get_wtime();
  pi = 12*(pi-SS);
  printf("Number of intervals = %d\n",n);
  printf(" Number of threads = %d\n", nthreads);
  printf(" pi is approximately %.16f\n", pi);
  printf(" Error is %e\n", fabs(pi-PI25DT));
  printf(" Wall clock time = %f\n", t1-t0);

  return 0;
}


struct Sd adapquad(double a0,double b0){
    double tol[10000];
    int n, m;
    double c, a[10000], b[10000], app[10000], Ints,oldapp;
    struct Sd s;
    Ints = 0.0;
    n = 0;
    m = 0;
    a[0] = a0;
    b[0] = b0;
    tol[0] = TOL;
    app[0] = trap(a[0],b[0]);
    m = m+1;
    while (n>=0){
      c = (a[n]+b[n])/2.0; oldapp = app[n];
      app[n] = trap(a[n],c); app[n+1] = trap(c,b[n]);
      if (fabs(oldapp-(app[n]+app[n+1]))<15*tol[n]){
         
         Ints = Ints +app[n]+app[n+1];
  /*       printf("%f\n",Ints);*/
         n = n-1;}
      else{
         tol[n] = tol[n]/2.0;
         tol[n+1] = tol[n];
         m = m+1;
         b[n+1] = b[n]; b[n] = c;
         a[n+1] = c;
         n = n+1;
      }
   } 
  s.value = 2*m+1;
  s.value2 = Ints;  
  return s;
} 


double trap(double a, double b){
double s, ao;
ao = (a+b)/2.0;
s = (sqrt(1.-a*a)+4*sqrt(1.-ao*ao)+sqrt(1.-b*b))*(b-a)/6.0;
  return s;
}
