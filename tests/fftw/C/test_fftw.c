#include <stdio.h>
#include <fftw3.h>

#define N 6
#define M 2

int main(int argc, char * argv[])
{
   double in[N*M] = {3.00, 1.50, -1.50, -3.00, -1.50, 1.50,
                     3.00, 1.50, -1.50, -3.00, -1.50, 1.50 }; // 3*cos(Th)
// double in[N*M] = {4.00, 0.00,  1.00,  0.00,  1.00, 0.00 }; // 1 + cos(1*Th) + cos(2*Th) + cos(3*Th)

//   fftw_complex out[(N/2+1)*(M+2)];
   fftw_complex * out;

   fftw_plan p_forward, p_reverse;
   int r, c, k;

   // Transform each column of a 2d array with 10 rows and 3 columns:
   //

   int rank    = 1;
   int howmany = 2;
   int istride = 1; /* distance between two elements in the same column */
   int ostride = 1;
   int idist   = N;
   int odist   = N/2+1;
   //int odist   = N;

   int n[] = {N};
   //   int inembed[] = {0};
   //   int onembed[] = {0};
   int * inembed = NULL;
   int * onembed = NULL;

   //in = (double*) fftw_malloc(sizeof(double) * N);
   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1) * M);

   printf("Creating plan...\n");

   //   p_forward = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);

   printf("Forward:\n");
   for (r = 0; r < M; r++) {
      for (c = 0; c < N; c++) {
         printf("%f ", in[c+r*idist]);
      }
      printf("\n");
   }
   printf("\n");

   p_forward = fftw_plan_many_dft_r2c(rank, n, howmany, 
                                      in , inembed, istride, idist,
                                      out, onembed, ostride, odist,
                                      FFTW_ESTIMATE);

   printf("Forward:\n");
   for (r = 0; r < M; r++) {
      for (c = 0; c < N; c++) {
         printf("%f ", in[c+r*idist]);
      }
      printf("\n");
   }
   printf("\n");

   fftw_execute_dft_r2c(p_forward, in, out);

   printf("Reverse:\n");
   for (r = 0; r < M; r++) {
      for (k = 0; k < N/2+1; k++) {
         printf("(%f %f) ", out[k+r*odist][0], out[k+r*odist][1]);
      }
      printf("\n");
   }
   printf("\n");

   for (k = 0; k < M*(N/2+1); k++) {
      printf("(%f %f) ", out[k][0], out[k][1]);
   }
   printf("\n");
   printf("\n");

   printf("Destroying plan... \n");

   fftw_destroy_plan(p_forward);
   //fftw_free(in);
   fftw_free(out);

#ifdef NOT_YET

   rank = 1;
 

   fftw_plan_many_dft_r2c

rfftwnd_plan rfftwnd_create_plan(int rank, const int *n,
                                 fftw_direction dir, int flags);

   p = rfftw_create_plan(N, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);

   rfftw_one(p, in, out);

   power_spectrum[0] = out[0]*out[0];  /* DC component */
   for (k = 1; k < (N+1)/2; ++k)  /* (k < N/2 rounded up) */
      power_spectrum[k] = out[k]*out[k] + out[N-k]*out[N-k];
   if (N % 2 == 0) /* N is even */
      power_spectrum[N/2] = out[N/2]*out[N/2];  /* Nyquist freq. */

      rfftw_destroy_plan(p);
#endif

   return 0;
}
