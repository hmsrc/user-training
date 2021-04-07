/*
 *  * OpenMP implementation of Monte Carlo pi-finding algorithm
 *   * (based on an example in Chap. 17 of Quinn (2004))
 *    *
 *     * usage: pi <samples> <threads>
 *      */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
int main (int argc, char *argv[])
{
  int i, count, samples, nthreads, seed;
  struct drand48_data drand_buf;
  double x, y;
  double t0, t1;

  samples  = atoi(argv[1]);
  nthreads = atoi(argv[2]);
  omp_set_num_threads (nthreads);

  t0 = omp_get_wtime();
  count = 0;

#pragma omp parallel private(i, x, y, seed, drand_buf) shared(samples)
{
  seed = 1202107158 + omp_get_thread_num() * 1999;
  srand48_r (seed, &drand_buf);

#pragma omp for reduction(+:count)
  for (i=0; i<samples; i++) {
    drand48_r (&drand_buf, &x);
    drand48_r (&drand_buf, &y);
    if (x*x + y*y <= 1.0) count++;
  }
}

  t1 = omp_get_wtime();
  printf("Estimate of pi: %7.5f\n", 4.0*count/samples);
  printf("Time: %7.2f\n", t1-t0);
}
