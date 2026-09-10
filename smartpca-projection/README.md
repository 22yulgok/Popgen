# smartpca persistent projection (EIGENSOFT v18140)

Extends EIGENSOFT **smartpca version 18140** so that a PCA computed once on a
reference dataset can be **saved to disk** and later **reused to project only
newly added individuals**, without recomputing the reference PCA.

The `.par` + `smartpca -p file.par` workflow is preserved. The feature adds a
small number of new `.par` parameters and one new module; it does **not** change
any of smartpca's PCA or projection mathematics.

## What this is (and is not)

- It reuses smartpca's **existing** code path. Projection of new samples goes
  through the **unmodified `lsqproj()`** routine, so genotype/SNP normalisation,
  allele-frequency handling, eigenvector scaling, reference-population handling,
  missing-genotype handling, and `lsqproject` behaviour are exactly stock.
- It does **not** use sklearn / PLINK / a generic `.evec` matrix multiply / any
  reimplementation or approximation of the projection.
- When none of the new parameters are set, smartpca behaves **exactly** as
  stock v18140 (verified byte-for-byte in the test suite).

## Repository layout

```
smartpca-projection/
├── README.md
├── environment.yml                 # conda build/test environment
├── scripts/
│   ├── build.sh                    # build the feature smartpca
│   └── run_tests.sh                # automated equivalence test
├── examples/
│   ├── reference.par               # Stage A: reference PCA + save model
│   ├── projection.par              # Stage B: load model + project new samples
│   └── run_pca_workflow.sh         # updated KKA-style two-stage shell workflow
├── tests/
│   ├── make_synth.py               # synthetic EIGENSTRAT data generator
│   └── compare_evec.py             # per-PC / max / mean .evec comparison
├── patch/
│   └── smartpca-persistent-projection.patch   # reviewable diff of smartpca.c
└── upstream/                       # EIGENSOFT source, pristine v18140 + feature
    ├── UPSTREAM_COMMIT.txt          # exact EIG commit the code is based on
    └── src/eigensrc/
        ├── smartpca.c               # v18140 + feature edits (annotated "ADDED")
        ├── smartpca.c.orig          # PRISTINE v18140 (for diffing / stock build)
        ├── smartpca_model.h         # NEW: model format + validation API
        └── smartpca_model.c         # NEW: model save/load + SNP alignment
```

`upstream/` is vendored from **DReichLab/EIG** at commit
`750e343b7f1800e54c1310c176c576e4ef52a136` (EIGENSOFT 8.0, `WVERSION "18140"`).
Upstream code is kept distinguishable from the modification in three ways:

1. `upstream/src/eigensrc/smartpca.c.orig` is the untouched v18140 file.
2. Every added region in `smartpca.c` is annotated with a
   `ADDED (feature/smartpca-persistent-projection)` comment.
3. `patch/smartpca-persistent-projection.patch` is the unified diff.

## Build

Inside the conda environment (recommended), the build uses the stock EIGENSOFT
dependencies (OpenBLAS/LAPACK + GSL), so the compiled binary matches upstream
numerics:

```bash
conda env create -f environment.yml
conda activate smartpca-projection

./scripts/build.sh          # -> build/smartpca
```

`build.sh` degrades gracefully:

- If GSL is unavailable, it disables only the fast-PCA (`fastmode`) path via a
  stub that aborts loudly if `fastmode` is ever requested. The classic
  `lsqproject` PCA and the projection feature are unaffected.
- If no system BLAS/LAPACK is found, it builds a static Reference-LAPACK
  locally so the build is self-contained (used for CI/sandbox).

## Test

```bash
./scripts/run_tests.sh
```

This builds both a pristine **stock** v18140 binary and the **feature** binary,
then runs the equivalence protocol below and prints a report.

### Equivalence protocol

**Baseline** — stock smartpca on `reference + test` individuals with
`poplistname` (reference pops) and `lsqproject: YES`; projected coordinates are
saved.

**Persistent-model test** —
1. Stage A: feature smartpca on the **reference only**, saving the model, then
   the process exits.
2. Stage B: a **new** feature-smartpca process loads only the saved model plus
   the test genotypes and projects them.

The Stage B coordinates are compared to the baseline for PCs 1–20, reporting the
**maximum**, **mean**, and **per-PC** absolute differences, for test individuals
with **complete**, **moderate (~30%)**, and **high (~80%)** missingness.

Because the saved basis is reused (not recomputed), eigenvector orientation is
preserved and no sign-flipping is applied.

### Observed result

On the synthetic data (3000 SNPs, 5 populations, 100 reference + 12 test):

- Backward compatibility: feature binary output is **byte-identical** (sha256)
  to stock when no model options are set.
- Projection equivalence: **max abs diff = 0.000e+00** across PCs 1–20 at a
  1e-6 tolerance (`hiprec`), including complete / moderate / high missingness.
- Allele-flip robustness: swapped-allele input is detected, recoded, and still
  matches to **0.000e+00**.
- Validation: strict mode (default) refuses to project when model SNPs are
  absent; `modelstrict: NO` projects from the available SNPs.

## Usage

### Stage A — compute reference PCA once and save the model

```
genotypename: reference.geno
snpname:      reference.snp
indivname:    reference.ind
poplistname:  reference.pops

evecoutname:  reference.evec
evaloutname:  reference.eval

altnormstype: NO
numoutevec:   20
numoutlieriter: 0
numoutlierevec: 0
outliersigmathresh: 6.0
numthreads:   8
qtmode:       0
lsqproject:   YES

modeloutname: reference.smartpca.model     # <-- NEW
```

```bash
smartpca -p reference.par
```

### Stage B — load the model and project new samples only

```
genotypename:      new_samples.geno
snpname:           new_samples.snp
indivname:         new_samples.ind

modelname:         reference.smartpca.model   # <-- NEW
projectionoutname: new_samples.evec           # <-- NEW
projectiononly:    YES                        # <-- NEW
# modelstrict:     YES  (default; NO allows model SNPs absent from new data)
```

```bash
smartpca -p projection.par
```

See `examples/run_pca_workflow.sh` for a two-stage version of the KKA workflow
that keeps the original loop/`sbatch` skeleton.

## New `.par` parameters

| Parameter            | Stage | Meaning |
|----------------------|-------|---------|
| `modeloutname`       | A     | Path to write the PCA projection model. Requires `lsqproject: YES`. |
| `modelname`          | B     | Path to a saved model to load. |
| `projectiononly`     | B     | `YES` to load the model and project new samples, then exit. |
| `projectionoutname`  | B     | `.evec` output path for the projected samples. |
| `modelstrict`        | B     | `YES` (default): every model SNP must be present in the new dataset. `NO`: absent model SNPs are treated as missing genotypes. |

## What the model stores

Determined from the actual smartpca source — only the state the existing
projection code consumes:

- SNP identity and ordering, allele coding (for validation + orientation)
- per-SNP normalisation: `xmean` (scaled mean) and `xfancy` (`1/sqrt(p(1-p))`)
- per-SNP × per-PC loadings `ffvecs`
- per-PC scaling constants `fxscal` and `eigscale`
- eigenvalues `lambda`, number of PCs `numeigs`
- normalisation configuration: `fancynorm`, `altnormstyle`, `missingmode`
- format magic + version and the smartpca `WVERSION` used

Doubles are written with full round-trip precision (`%.17g`).

### Model validation

`read_pca_model` + `model_check_and_align` validate compatibility before
projecting, and **fail explicitly** when correctness cannot be guaranteed:

- model format magic + version
- number of PCs and normalisation configuration
- SNP identity and correspondence (matched by SNP ID)
- allele orientation (identical, swap→recode `g→2-g`, or a fatal conflict)

SNP rows are **never** silently assumed to correspond.

## Notes / limitations

- The equivalence tests use synthetic EIGENSTRAT data so they run without any
  private data or server paths.
- The `fastmode` (randomized fast PCA) code path requires GSL and is orthogonal
  to this feature; the sandbox/CI build stubs it out.
