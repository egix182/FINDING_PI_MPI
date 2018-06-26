#include <stdio.h>
#include <time.h>

#define N 1E7
#define d 1E-7
#define d2 1E-14

int main (int argc, char* argv[])
{
   clock_t begin = clock(), end;
   double pi=0.0, result=0.0,x2=0.0;
   int i=0;

   for (i=0; i<N; i+=1)
   {
       x2=d2*i*i;
       result+=1.0/(1.0+x2);
   }
   
   pi=4*d*result;
   end = clock(); 
   printf("PI = %lf\t Time = %lf\n", pi, ((double)(end-begin))/CLOCKS_PER_SEC);
   
   return 0;
}
