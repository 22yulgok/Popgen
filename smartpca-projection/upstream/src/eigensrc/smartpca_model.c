/*
 * smartpca_model.c  --  persistent PCA projection model I/O + validation.
 *
 * ADDED for the "persistent projection" feature.  NEW code; does not change
 * any stock PCA/projection mathematics.  See smartpca_model.h for the list of
 * persisted quantities and the rationale.
 *
 * On-disk format (plain text, self-describing, auditable):
 *
 *   SMARTPCA-PROJMODEL
 *   format_version <int>
 *   smartpca_version <string>
 *   ncols <int>
 *   numeigs <int>
 *   nrows <int>
 *   fancynorm <int>
 *   altnormstyle <int>
 *   missingmode <int>
 *   lambda: v0 v1 ... v(numeigs-1)
 *   fxscal: v0 ... v(numeigs-1)
 *   eigscale: v0 ... v(numeigs-1)
 *   SNPS <ncols>
 *   <id> <chrom> <physpos> <a0> <a1> <xmean> <xfancy> <load0> <load1> ... <load(numeigs-1)>
 *   ... (ncols such lines) ...
 *   END
 *
 * All doubles are written with %.17g so that IEEE-754 doubles round-trip
 * exactly; the reloaded model therefore reproduces stock smartpca numbers to
 * the bit, not merely to a tolerance.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strsubs.h>            /* fatalx, openit, ZALLOC */

#include "smartpca_model.h"

/* ---------------------------------------------------------------- writing */

int
write_pca_model (const char *path,
                 const char *smartpca_version,
                 int ncols, int numeigs, int nrows,
                 int fancynorm, int altnormstyle, int missingmode,
                 SNP ** xsnplist,
                 double *xmean, double *xfancy,
                 double *ffvecs, double *fxscal, double *eigscale,
                 double *lambda)
{
  FILE *fp;
  int i, j;
  SNP *cupt;

  openit ((char *) path, &fp, "w");     /* fatalx on failure */

  fprintf (fp, "%s\n", SMARTPCA_MODEL_MAGIC);
  fprintf (fp, "format_version %d\n", SMARTPCA_MODEL_VERSION);
  fprintf (fp, "smartpca_version %s\n", smartpca_version);
  fprintf (fp, "ncols %d\n", ncols);
  fprintf (fp, "numeigs %d\n", numeigs);
  fprintf (fp, "nrows %d\n", nrows);
  fprintf (fp, "fancynorm %d\n", fancynorm);
  fprintf (fp, "altnormstyle %d\n", altnormstyle);
  fprintf (fp, "missingmode %d\n", missingmode);

  fprintf (fp, "lambda:");
  for (j = 0; j < numeigs; ++j)
    fprintf (fp, " %.17g", lambda[j]);
  fprintf (fp, "\n");

  fprintf (fp, "fxscal:");
  for (j = 0; j < numeigs; ++j)
    fprintf (fp, " %.17g", fxscal[j]);
  fprintf (fp, "\n");

  fprintf (fp, "eigscale:");
  for (j = 0; j < numeigs; ++j)
    fprintf (fp, " %.17g", eigscale[j]);
  fprintf (fp, "\n");

  fprintf (fp, "SNPS %d\n", ncols);
  for (i = 0; i < ncols; ++i) {
    cupt = xsnplist[i];
    fprintf (fp, "%s %d %.17g %c %c %.17g %.17g",
             cupt->ID, cupt->chrom, cupt->physpos,
             cupt->alleles[0] ? cupt->alleles[0] : '?',
             cupt->alleles[1] ? cupt->alleles[1] : '?',
             xmean[i], xfancy[i]);
    for (j = 0; j < numeigs; ++j)
      fprintf (fp, " %.17g", ffvecs[j * ncols + i]);
    fprintf (fp, "\n");
  }
  fprintf (fp, "END\n");

  fclose (fp);
  printf ("## wrote PCA projection model: %s  (ncols=%d numeigs=%d nrows=%d)\n",
          path, ncols, numeigs, nrows);
  return 0;
}

/* ---------------------------------------------------------------- reading */

static void
expect_key (FILE * fp, const char *key)
{
  char buf[256];
  if (fscanf (fp, "%255s", buf) != 1)
    fatalx ("model file: unexpected EOF, wanted '%s'\n", key);
  if (strcmp (buf, key) != 0)
    fatalx ("model file: expected key '%s' but found '%s'\n", key, buf);
}

PcaModel *
read_pca_model (const char *path)
{
  FILE *fp;
  PcaModel *m;
  char buf[256];
  int i, j;

  openit ((char *) path, &fp, "r");     /* fatalx on failure */

  if (fscanf (fp, "%255s", buf) != 1)
    fatalx ("model file: empty file %s\n", path);
  if (strcmp (buf, SMARTPCA_MODEL_MAGIC) != 0)
    fatalx ("model file: bad magic '%s' (not a smartpca projection model)\n",
            buf);

  ZALLOC (m, 1, PcaModel);

  expect_key (fp, "format_version");
  if (fscanf (fp, "%d", &m->format_version) != 1)
    fatalx ("model file: bad format_version\n");
  if (m->format_version != SMARTPCA_MODEL_VERSION)
    fatalx ("model file: format_version %d != supported %d\n",
            m->format_version, SMARTPCA_MODEL_VERSION);

  expect_key (fp, "smartpca_version");
  if (fscanf (fp, "%15s", m->smartpca_version) != 1)
    fatalx ("model file: bad smartpca_version\n");

  expect_key (fp, "ncols");
  if (fscanf (fp, "%d", &m->ncols) != 1)
    fatalx ("model file: bad ncols\n");
  expect_key (fp, "numeigs");
  if (fscanf (fp, "%d", &m->numeigs) != 1)
    fatalx ("model file: bad numeigs\n");
  expect_key (fp, "nrows");
  if (fscanf (fp, "%d", &m->nrows) != 1)
    fatalx ("model file: bad nrows\n");
  expect_key (fp, "fancynorm");
  if (fscanf (fp, "%d", &m->fancynorm) != 1)
    fatalx ("model file: bad fancynorm\n");
  expect_key (fp, "altnormstyle");
  if (fscanf (fp, "%d", &m->altnormstyle) != 1)
    fatalx ("model file: bad altnormstyle\n");
  expect_key (fp, "missingmode");
  if (fscanf (fp, "%d", &m->missingmode) != 1)
    fatalx ("model file: bad missingmode\n");

  if (m->ncols <= 0 || m->numeigs <= 0)
    fatalx ("model file: nonsensical dimensions ncols=%d numeigs=%d\n",
            m->ncols, m->numeigs);

  ZALLOC (m->lambda, m->numeigs, double);
  ZALLOC (m->fxscal, m->numeigs, double);
  ZALLOC (m->eigscale, m->numeigs, double);

  expect_key (fp, "lambda:");
  for (j = 0; j < m->numeigs; ++j)
    if (fscanf (fp, "%lf", &m->lambda[j]) != 1)
      fatalx ("model file: bad lambda[%d]\n", j);
  expect_key (fp, "fxscal:");
  for (j = 0; j < m->numeigs; ++j)
    if (fscanf (fp, "%lf", &m->fxscal[j]) != 1)
      fatalx ("model file: bad fxscal[%d]\n", j);
  expect_key (fp, "eigscale:");
  for (j = 0; j < m->numeigs; ++j)
    if (fscanf (fp, "%lf", &m->eigscale[j]) != 1)
      fatalx ("model file: bad eigscale[%d]\n", j);

  expect_key (fp, "SNPS");
  {
    int nsnp;
    if (fscanf (fp, "%d", &nsnp) != 1)
      fatalx ("model file: bad SNP count\n");
    if (nsnp != m->ncols)
      fatalx ("model file: SNP count %d != ncols %d\n", nsnp, m->ncols);
  }

  ZALLOC (m->snps, m->ncols, ModelSNP);
  ZALLOC (m->xmean, m->ncols, double);
  ZALLOC (m->xfancy, m->ncols, double);
  ZALLOC (m->ffvecs, (long) m->ncols * m->numeigs, double);

  for (i = 0; i < m->ncols; ++i) {
    char a0[8], a1[8];
    ModelSNP *s = &m->snps[i];
    if (fscanf (fp, "%39s %d %lf %7s %7s %lf %lf",
                s->snpid, &s->chrom, &s->physpos, a0, a1,
                &m->xmean[i], &m->xfancy[i]) != 7)
      fatalx ("model file: bad SNP record at line %d\n", i);
    s->alleles[0] = a0[0];
    s->alleles[1] = a1[0];
    for (j = 0; j < m->numeigs; ++j)
      if (fscanf (fp, "%lf", &m->ffvecs[j * m->ncols + i]) != 1)
        fatalx ("model file: bad loading SNP %d PC %d\n", i, j);
  }

  expect_key (fp, "END");
  fclose (fp);

  printf ("## loaded PCA projection model: %s\n", path);
  printf ("##   smartpca_version=%s ncols=%d numeigs=%d nrows=%d\n",
          m->smartpca_version, m->ncols, m->numeigs, m->nrows);
  printf ("##   fancynorm=%d altnormstyle=%d missingmode=%d\n",
          m->fancynorm, m->altnormstyle, m->missingmode);
  return m;
}

void
free_pca_model (PcaModel * m)
{
  if (m == NULL)
    return;
  if (m->snps) free (m->snps);
  if (m->xmean) free (m->xmean);
  if (m->xfancy) free (m->xfancy);
  if (m->ffvecs) free (m->ffvecs);
  if (m->fxscal) free (m->fxscal);
  if (m->eigscale) free (m->eigscale);
  if (m->lambda) free (m->lambda);
  free (m);
}

/* ---------------------------------------------------- validate + align */

static int
allele_ok (char m0, char m1, char n0, char n1, int *flip)
{
  /* Unknown alleles ('?', '0', '1'/'2' placeholder coding): accept without
   * flipping (dosage coding is assumed identical, as in stock EIGENSTRAT). */
  *flip = 0;
  if (m0 == '?' || m1 == '?' || n0 == '?' || n1 == '?')
    return 1;
  /* EIGENSTRAT numeric placeholder coding '1'/'2' carries no strand info. */
  if ((m0 == '1' && m1 == '2') || (n0 == '1' && n1 == '2'))
    return 1;
  if (m0 == n0 && m1 == n1) {
    *flip = 0;
    return 1;
  }
  if (m0 == n1 && m1 == n0) {
    *flip = 1;                  /* alleles swapped: recode dosage g -> 2-g */
    return 1;
  }
  return 0;                     /* genuinely incompatible alleles */
}

int
model_check_and_align (PcaModel * m,
                       SNP ** newsnps, int numnewsnps,
                       int strict,
                       int **out_order, int **out_flip)
{
  int *order, *flip;
  int i, k;
  int nfound = 0, nflip = 0, nmiss = 0, nbadallele = 0;
  /* Build a lookup from SNP id -> index in the new dataset. */
  /* Linear scan is O(ncols*numnewsnps); fine for validation and keeps the
   * code dependency-free.  A hash could replace this if needed. */

  ZALLOC (order, m->ncols, int);
  ZALLOC (flip, m->ncols, int);
  for (i = 0; i < m->ncols; ++i)
    order[i] = -1;

  for (i = 0; i < m->ncols; ++i) {
    ModelSNP *s = &m->snps[i];
    int hit = -1;
    for (k = 0; k < numnewsnps; ++k) {
      if (strcmp (newsnps[k]->ID, s->snpid) == 0) {
        hit = k;
        break;
      }
    }
    if (hit < 0) {
      ++nmiss;
      if (strict) {
        if (nmiss <= 10)
          printf ("## MODEL ERROR: SNP '%s' (model col %d) absent from new dataset\n",
                  s->snpid, i);
        continue;
      }
      order[i] = -1;            /* treated as missing for every individual */
      continue;
    }
    {
      int fl = 0;
      SNP *ns = newsnps[hit];
      if (!allele_ok (s->alleles[0], s->alleles[1],
                      ns->alleles[0], ns->alleles[1], &fl)) {
        ++nbadallele;
        if (nbadallele <= 10)
          printf ("## MODEL ERROR: SNP '%s' allele mismatch model=%c/%c new=%c/%c\n",
                  s->snpid, s->alleles[0], s->alleles[1],
                  ns->alleles[0], ns->alleles[1]);
        continue;
      }
      order[i] = hit;
      flip[i] = fl;
      ++nfound;
      if (fl)
        ++nflip;
    }
  }

  printf ("## model/dataset SNP check: matched=%d flipped=%d missing=%d allele_conflicts=%d (model ncols=%d)\n",
          nfound, nflip, nmiss, nbadallele, m->ncols);

  if (nbadallele > 0) {
    free (order);
    free (flip);
    fatalx ("model/dataset allele conflicts (%d): refusing to project (correctness cannot be guaranteed)\n",
            nbadallele);
  }
  if (strict && nmiss > 0) {
    free (order);
    free (flip);
    fatalx ("model/dataset SNP mismatch: %d model SNPs absent from new dataset (strict mode)\n",
            nmiss);
  }
  if (nfound == 0) {
    free (order);
    free (flip);
    fatalx ("model/dataset SNP check: no SNPs matched by ID; cannot project\n");
  }

  *out_order = order;
  *out_flip = flip;
  return 0;
}
