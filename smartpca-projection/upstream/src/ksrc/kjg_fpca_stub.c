/*
 * kjg_fpca_stub.c  --  TEST-BUILD-ONLY stub for the fast randomized PCA path.
 *
 * The stock smartpca fast-PCA path (fastmode: YES) lives in ksrc/kjg_fpca.c and
 * depends on GSL.  The persistent-projection feature and its equivalence tests
 * use ONLY the classic (lsqproject) PCA path, which never calls kjg_fpca().
 *
 * This stub lets the sandbox build link without GSL.  It is NOT part of the
 * algorithm and is NOT used by any tested code path.  A production build in the
 * conda environment links the real ksrc/kjg_fpca.o against GSL instead.
 *
 * If fastmode is ever requested against a stub-linked binary, we abort loudly
 * rather than silently returning wrong numbers.
 */
#include <stdio.h>
#include <stdlib.h>

void
kjg_fpca (int K, int L, int I, double *evals, double *evecs)
{
  (void) K; (void) L; (void) I; (void) evals; (void) evecs;
  fprintf (stderr,
           "FATAL: this smartpca binary was built WITHOUT the GSL fast-PCA "
           "path (fastmode). Rebuild against GSL (ksrc/kjg_fpca.o) to use "
           "fastmode. The classic lsqproject PCA path does not need this.\n");
  exit (1);
}
