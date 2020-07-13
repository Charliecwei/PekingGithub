/* 
  Foundations of Parallel and Distributed Computing,Homework 3 Falls 2019.
  Instructor: Chen Wei@ Peking University.
  This code shows how to calculate areas of Mandelbrot Set.
*/


#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>




#define Y 1.5065918849  //Areas of Mandelbrot Set. 
//#define X_Grid_N 1000    //Mesh fraction
//#define Y_Grid_N 1000   //Mesh fraction
//#define Max_iter 5000 //Maximum number of iterations
struct complex{double x;double y;};





static int Mandelbrot_Setnode(double a, double b, int Max_iter){
	int i, ns;
	struct complex c, z, zold;
	z.x = 0;
	z.y = 0;
	zold = z;
	c.x = a;
	c.y = b;
	//k<Max_iter；
	for (i = 0; i < Max_iter; i++){
		z.x = zold.x * zold.x - zold.y * zold.y + c.x;
		z.y = 2.0 * zold.x * zold.y + c.y;
		zold = z;
		//printf(" z.x = %f, z.y = %f\n", z.x, z.y);
		if  ((z.x * z.x + z.y * z.y) >= 4.0)             //(abs(z)>=2)
            break;
	}

	 int s = 0;
	 if (i == Max_iter){
	 	s = 1;
	 }
	 return s;
}






static double Mandelbrot_SetAreas(double a, double b, double hx, double hy){
 double Area = 0; //The area of the Mandelbrit set cube
// srand((unsigned)time(NULL));
 srand((unsigned)1);
 int i, ns, n = 0, m = 0, Max_iter = 1536, s = 0;  //1549,1525,1530,1535,
 double lax,lay, hxmin, hymin;
 double ax[1000], bx[1000],ay[1000],by[1000];
 ax[0] = a;
 bx[0] = a+hx;
 ay[0] = b;
 by[0] = b+hy;
 hymin = hy/32.0;
 hxmin = hx/32.0;
 
while (n>=0){
	for (i=0; i<10; i++){
				//lax = rand()/(RAND_MAX+1.0);
				//lay = rand()/(RAND_MAX+1.0);
		         lax = ((double)i+0.5)/12.;
		         lay = ((double)i+0.5)/12.0;
				s += Mandelbrot_Setnode(ax[n]+lax*hx, ay[n]+lay*hy, Max_iter);
			}
		if (s>=7||s<=3||(hy <= hymin)||(hx<=hxmin)){
		Area = Area + hx*hy*(double)s/10.0;
		n = n-1;
		Max_iter = Max_iter/2;
		
	    }

		else {
			double cx = (ax[n]+bx[n])/2.0;
			double cy = (ay[n]+by[n])/2.0;
		    ax[n+1] = ax[n];
		    ax[n+2] = cx;
		    ax[n+3] = cx;
		    bx[n+3] = bx[n];
		    bx[n+2] = bx[n];
		    bx[n+1] = cx;
		    bx[n] = cx;

		    ay[n+1] = cy;
		    ay[n+2] = ay[n];
		    ay[n+3] = cy;
		    by[n+3] = by[n];
		    by[n+2] = cy;
		    by[n+1] = by[n];
		    by[n] = cy;
		    n = n + 1;
		 
		    	  
        }
       // printf("%d\n", n);
        hx = (bx[n] - ax[n]);
		hy = (by[n] - ay[n]);
		s = 0;
		if (Max_iter<50000){
		Max_iter = 2*Max_iter;
       }
}



return Area;
	
}






int main(int argc, char *argv[]){

 FILE *fpWrite=fopen("data.txt","w");//wirite data

if(fpWrite==NULL)

{return 0;

}




    
	
	double a, b, t0, t1;
	int i,j,ii,jj, nthreads;
	int X_Grid_N, Y_Grid_N;
	int X_GRIDE[5] = {800,1000,1200,1400,1600};
	int Y_GRIDE[5] = {800,1000,1200,1400,1600};

for (ii = 1; ii<2; ii++){
	for (jj = 1; jj<2; jj++){
		X_Grid_N = X_GRIDE[ii];
		Y_Grid_N = Y_GRIDE[jj];
		
		double Areas = 0;  //Area approximation of Mandelbrot Set.
		double hx = 3.0/(double) X_Grid_N; //The grid size of x, x\in [-2.0,0.8];
		double hy = 1.5/(double) Y_Grid_N; //The grid size of y, y\in [0,1.5];   Half of the Mandelbrot Set area;
   
			t0 = omp_get_wtime();
			nthreads = omp_get_max_threads();
			#pragma omp parallel for private(i, j, a, b) reduction(+:Areas)  schedule(dynamic) 
			for (i=0; i<=X_Grid_N-1; i++){
				a = -2.0+(double) i * hx;
					for (j=0; j<=Y_Grid_N-1; j++){
			            b = (double) j * hy;
			         
			            Areas += Mandelbrot_SetAreas(a, b, hx, hy); //Adaptive algorithm is used to calculate area of [a,a+hx]\times [b,b+hy];

					}
			}


				Areas = 2.0 * Areas;
				t1  = omp_get_wtime();
			
			printf(" Number of threads = %d\n", nthreads);
			printf(" Mandelbrot Areas is %.9f\n", Y);	
		    printf(" Areas is approximately  is %.9f\n", Areas);
			printf(" Error is %.5e\n", fabs(Y - Areas));
			printf(" Wall clock time = %f\n", t1-t0);

		    fprintf(fpWrite,"%.5e\n ",fabs(Y - Areas));



	}
}
	
fclose(fpWrite);//	
return 0;

}
