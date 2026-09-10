/*
 * smartpca_model.h  --  persistent PCA projection model for smartpca.
 *
 * ADDED for the "persistent projection" feature (feature/smartpca-persistent-
 * projection).  This is NEW code layered on top of stock EIGENSOFT smartpca
 * v18140; it does not alter any PCA/projection mathematics.  It only
 * serialises / deserialises the internal state that stock smartpca already
 * computes (per-SNP normalisation, SNP loadings, per-PC scaling constants),
 * and validates that a saved model is compatible with a new dataset.
 *
 * The model captures exactly the quantities the stock lsqproject code path
 * consumes when projecting an individual:
 *
 *   xmean[ncols]            per-SNP scaled mean               (fixxrow)
 *   xfancy[ncols]           per-SNP normalisation multiplier  (fixxrow)
 *   ffvecs[ncols*numeigs]   per-SNP x per-PC loadings         (emat / dot)
 *   fxscal[numeigs]         per-PC dot-product scaling        (emat / dot)
 *   eigscale[numeigs]       per-PC lsq->dot rescale           (mulmat)
 *   lambda[numeigs]         eigenvalues (for the .evec header / .eval)
 *   per-SNP identity        ID, chrom, physpos, alleles[2]    (validation)
 *   normalisation flags     fancynorm, altnormstyle, missingmode
 *
 * With these, projecting a NEW individual is exactly the stock computation
 *   loadxdataind -> fixxrow(xmean,xfancy) -> regressit over observed SNPs
 *   -> multiply by eigscale
 * performed by the unmodified lsqproj() routine, with no recomputation of the
 * reference eigenanalysis.
 */

#ifndef SMARTPCA_MODEL_H
#define SMARTPCA_MODEL_H

#include "admutils.h"

/* Bump only on incompatible on-disk format changes. */
#define SMARTPCA_MODEL_MAGIC   "SMARTPCA-PROJMODEL"
#define SMARTPCA_MODEL_VERSION 1

typedef struct
{
  char snpid[IDSIZE];
  int chrom;
  double physpos;
  char alleles[2];              /* reference/variant allele coding */
} ModelSNP;

typedef struct
{
  /* format / provenance */
  int format_version;           /* SMARTPCA_MODEL_VERSION */
  char smartpca_version[16];    /* WVERSION the model was written with */

  /* dimensions */
  int ncols;                    /* number of SNPs in the reference PCA */
  int numeigs;                  /* number of PCs */
  int nrows;                    /* number of reference individuals (info) */

  /* normalisation configuration (must match to reproduce numbers) */
  int fancynorm;
  int altnormstyle;
  int missingmode;

  /* per-SNP state (length ncols) */
  ModelSNP *snps;
  double *xmean;
  double *xfancy;

  /* projection basis */
  double *ffvecs;               /* ncols*numeigs, PC-major: ffvecs[j*ncols+k] */
  double *fxscal;               /* numeigs */
  double *eigscale;             /* numeigs */
  double *lambda;               /* numeigs (top eigenvalues) */
} PcaModel;

/*
 * Serialise a model.  Writes a small self-describing text header followed by
 * the SNP table and the numeric arrays.  Returns 0 on success.
 *
 * Called by smartpca once the reference PCA is complete (ffvecs, fxscal,
 * eigscale, lambda all populated).
 */
int write_pca_model (const char *path,
                     const char *smartpca_version,
                     int ncols, int numeigs, int nrows,
                     int fancynorm, int altnormstyle, int missingmode,
                     SNP ** xsnplist,
                     double *xmean, double *xfancy,
                     double *ffvecs, double *fxscal, double *eigscale,
                     double *lambda);

/* Load a model from disk.  Returns a malloc'd PcaModel* (free with
 * free_pca_model) or NULL on failure. */
PcaModel *read_pca_model (const char *path);

void free_pca_model (PcaModel * m);

/*
 * Validate that a loaded model is compatible with the new dataset's SNP set,
 * and (as a side-effect) reorder / mask the new dataset so its SNP columns
 * correspond exactly to the model's SNP ordering.
 *
 * newsnps/numnewsnps describe the SNPs loaded from the projection dataset.
 * On success, *out_order is a malloc'd int[ncols] giving, for each model SNP
 * column, the index into newsnps[] to use (or -1 if that SNP is absent from
 * the new dataset and must be treated as fully missing).  *out_flip is a
 * malloc'd int[ncols]: 1 if the new dataset's alleles are swapped relative to
 * the model (dosage must be recoded g -> 2-g), else 0.
 *
 * Fails explicitly (returns non-zero, prints reason) when correctness cannot
 * be guaranteed: SNP identity mismatch beyond tolerance, allele incompatibility,
 * numeigs/normalisation/version mismatch.  strict controls whether a SNP
 * missing from the new dataset is fatal (strict=1) or allowed-as-missing
 * (strict=0).
 */
int model_check_and_align (PcaModel * m,
                           SNP ** newsnps, int numnewsnps,
                           int strict,
                           int **out_order, int **out_flip);

#endif /* SMARTPCA_MODEL_H */
