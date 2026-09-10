#!/usr/bin/env bash
#
# run_tests.sh -- automated equivalence test for the persistent-projection
# feature.
#
# It builds TWO binaries:
#   * smartpca_stock  -- pristine EIGENSOFT v18140 (from smartpca.c.orig), used
#                        to produce the BASELINE.
#   * smartpca        -- the feature build (from build.sh).
#
# Then:
#   Baseline        : stock smartpca on (reference + test) with poplistname +
#                     lsqproject: YES  -> baseline.evec (projected coords).
#   Stage A         : feature smartpca on (reference only), saves a model.
#   [process exits]
#   Stage B         : a NEW feature-smartpca process loads ONLY the model + the
#                     test genotypes and projects them -> stageB.evec.
#
# Finally it compares Stage B against the baseline for PCs 1..20 and reports the
# maximum / mean / per-PC absolute differences, for test individuals with
# complete, moderate, and high genotype missingness.  It also checks backward
# compatibility (feature binary == stock binary when no model options are set).
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRCROOT="$HERE/upstream"
SRC="$SRCROOT/src"
BUILD="$HERE/build"
WORK="${WORK:-$HERE/build/testrun}"
TESTS="$HERE/tests"
TOL="${TOL:-1e-6}"
NPC="${NPC:-20}"

CC="${CC:-cc}"
PY="${PYTHON:-python3}"

mkdir -p "$WORK"

# ---------------------------------------------------------------- 1. build
echo "== [1/5] build feature smartpca =="
"$HERE/scripts/build.sh"

echo "== [2/5] build pristine stock smartpca (baseline) =="
# Compile the pristine smartpca.c.orig into a stock binary, reusing the objects
# and libraries already produced by build.sh.
CFLAGS="-O2 -g -Wimplicit -I$SRCROOT/include -I$SRC/nicksrc"
$CC $CFLAGS -x c -c -o "$WORK/smartpca_stock.o" "$SRC/eigensrc/smartpca.c.orig"

# Reuse the same link line the feature build used, minus the model object.
# Detect libs the same way build.sh did (BLAS + optional GSL).
BLAS_LIBS=""; GSL_LIBS=""; KOBJS="$SRC/ksrc/kjg_fpca_stub.o"
if [ -f "$BUILD/lapack/liblapack.a" ]; then
  BLAS_LIBS="$BUILD/lapack/liblapack.a $BUILD/lapack/librefblas.a -lgfortran"
elif [ -n "${CONDA_PREFIX:-}" ] && ls "$CONDA_PREFIX"/lib/libopenblas* >/dev/null 2>&1; then
  BLAS_LIBS="-L$CONDA_PREFIX/lib -lopenblas"
elif ls /usr/lib*/libopenblas* /usr/lib/*/libopenblas* >/dev/null 2>&1; then
  BLAS_LIBS="-lopenblas"
else
  BLAS_LIBS="-llapack -lblas"
fi
if [ -f "$SRC/ksrc/kjg_fpca.o" ]; then
  KOBJS="$SRC/ksrc/kjg_fpca.o $SRC/ksrc/kjg_gsl.o"
  command -v pkg-config >/dev/null 2>&1 && pkg-config --exists gsl 2>/dev/null && GSL_LIBS="$(pkg-config --libs gsl)"
fi

$CC -O2 -g -o "$WORK/smartpca_stock" \
  "$WORK/smartpca_stock.o" \
  "$SRC/eigensrc/eigsubs.o" "$SRC/eigensrc/exclude.o" "$SRC/eigensrc/smartsubs.o" "$SRC/eigensrc/eigx.o" \
  "$SRC/mcio.o" "$SRC/qpsubs.o" "$SRC/admutils.o" "$SRC/egsubs.o" "$SRC/regsubs.o" "$SRC/gval.o" \
  $KOBJS "$SRC/nicksrc/libnick.a" \
  $BLAS_LIBS $GSL_LIBS -lm -lpthread

STOCK="$WORK/smartpca_stock"
FEAT="$BUILD/smartpca"
"$STOCK" -v 2>&1 | grep -i version | head -1

# ---------------------------------------------------------------- 2. data
echo "== [3/5] generate synthetic EIGENSTRAT data =="
DATA="$WORK/data"
"$PY" "$TESTS/make_synth.py" --out "$DATA" \
  --nsnps "${NSNPS:-3000}" --npops "${NPOPS:-5}" \
  --per-pop "${PERPOP:-20}" --ntest "${NTEST:-12}" --seed "${SEED:-42}"

gen_par () {  # helper: emit a .par file
  local f="$1"; shift; : > "$f"; for kv in "$@"; do echo "$kv" >> "$f"; done
}

# ---------------------------------------------------------------- 3. baseline
echo "== [4/5] baseline (stock, reference+test, lsqproject) =="
gen_par "$WORK/baseline.par" \
  "genotypename: $DATA/all.geno" \
  "snpname: $DATA/all.snp" \
  "indivname: $DATA/all.ind" \
  "poplistname: $DATA/ref.pops" \
  "evecoutname: $WORK/baseline.evec" \
  "evaloutname: $WORK/baseline.eval" \
  "altnormstype: NO" \
  "numoutevec: $NPC" \
  "numoutlieriter: 0" \
  "numoutlierevec: 0" \
  "outliersigmathresh: 6.0" \
  "numthreads: 8" \
  "qtmode: 0" \
  "lsqproject: YES" \
  "hiprec: YES"
"$STOCK" -p "$WORK/baseline.par" > "$WORK/baseline.log" 2>&1

# Backward-compatibility: feature binary must produce identical output with no
# model options set.
"$FEAT" -p "$WORK/baseline.par" > "$WORK/baseline_feat.log" 2>&1
cp "$WORK/baseline.evec" "$WORK/baseline_feat.evec"
"$STOCK" -p "$WORK/baseline.par" > /dev/null 2>&1   # restore stock output

# ---------------------------------------------------------------- 4. A + B
echo "== Stage A (feature, reference only, save model) =="
gen_par "$WORK/stageA.par" \
  "genotypename: $DATA/ref.geno" \
  "snpname: $DATA/ref.snp" \
  "indivname: $DATA/ref.ind" \
  "poplistname: $DATA/ref.pops" \
  "evecoutname: $WORK/stageA.evec" \
  "evaloutname: $WORK/stageA.eval" \
  "altnormstype: NO" \
  "numoutevec: $NPC" \
  "numoutlieriter: 0" \
  "numoutlierevec: 0" \
  "outliersigmathresh: 6.0" \
  "numthreads: 8" \
  "qtmode: 0" \
  "lsqproject: YES" \
  "modeloutname: $WORK/reference.smartpca.model"
"$FEAT" -p "$WORK/stageA.par" > "$WORK/stageA.log" 2>&1

echo "== Stage B (NEW process, load model, project test only) =="
gen_par "$WORK/stageB.par" \
  "genotypename: $DATA/test.geno" \
  "snpname: $DATA/test.snp" \
  "indivname: $DATA/test.ind" \
  "modelname: $WORK/reference.smartpca.model" \
  "projectionoutname: $WORK/stageB.evec" \
  "projectiononly: YES" \
  "hiprec: YES"
"$FEAT" -p "$WORK/stageB.par" > "$WORK/stageB.log" 2>&1

# ---------------------------------------------------------------- 5. compare
echo "== [5/5] equivalence report =="
rc=0

echo
echo "---- backward compatibility (feature vs stock, no model options) ----"
if "$PY" - "$WORK/baseline.evec" "$WORK/baseline_feat.evec" <<'EOF'
import sys,hashlib
a,b=sys.argv[1],sys.argv[2]
ha=hashlib.sha256(open(a,'rb').read()).hexdigest()
hb=hashlib.sha256(open(b,'rb').read()).hexdigest()
print("stock .evec  sha256:",ha[:16])
print("feature .evec sha256:",hb[:16])
sys.exit(0 if ha==hb else 1)
EOF
then echo "backward compatibility: PASS (byte-identical)"; else echo "backward compatibility: FAIL"; rc=1; fi

echo
echo "---- projection equivalence (Stage B vs stock baseline), all samples ----"
"$PY" "$TESTS/compare_evec.py" "$WORK/baseline.evec" "$WORK/stageB.evec" --tol "$TOL" --npc "$NPC" || rc=1

# The synthetic test IDs are tagged COMPLETE / MOD## / HIGH## by missingness.
for tier in COMPLETE MOD HIGH; do
  echo
  echo "---- missingness tier: $tier ----"
  grep -h "TEST_.*_${tier}" "$WORK/baseline.evec" > "$WORK/base_${tier}.evec" || true
  grep -h "TEST_.*_${tier}" "$WORK/stageB.evec"   > "$WORK/proj_${tier}.evec" || true
  if [ -s "$WORK/base_${tier}.evec" ]; then
    "$PY" "$TESTS/compare_evec.py" "$WORK/base_${tier}.evec" "$WORK/proj_${tier}.evec" \
      --tol "$TOL" --npc "$NPC" | grep -E "compared|maximum abs|mean abs|RESULT" || rc=1
    if "$PY" "$TESTS/compare_evec.py" "$WORK/base_${tier}.evec" "$WORK/proj_${tier}.evec" \
        --tol "$TOL" --npc "$NPC" >/dev/null 2>&1; then :; else rc=1; fi
  fi
done

# ---- allele-orientation robustness: swap alleles + recode dosages, project,
#      and confirm the result is unchanged.
echo
echo "---- allele-flip robustness (swapped alleles must project identically) ----"
"$PY" - "$DATA/test.snp" "$DATA/test.geno" "$DATA/test.ind" "$WORK/flip" <<'EOF'
import sys
snp,geno,ind,out=sys.argv[1:5]
L=open(snp).read().splitlines()
o=[]
for line in L:
    p=line.split(); p[4],p[5]=p[5],p[4]
    o.append("%20s %2s %12s %12s %s %s"%(p[0],p[1],p[2],p[3],p[4],p[5]))
open(out+".snp","w").write("\n".join(o)+"\n")
rec={'0':'2','1':'1','2':'0','9':'9'}
g=open(geno).read().splitlines()
open(out+".geno","w").write("\n".join("".join(rec[c] for c in r) for r in g)+"\n")
open(out+".ind","w").write(open(ind).read())
EOF
gen_par "$WORK/flip.par" \
  "genotypename: $WORK/flip.geno" \
  "snpname: $WORK/flip.snp" \
  "indivname: $WORK/flip.ind" \
  "modelname: $WORK/reference.smartpca.model" \
  "projectionoutname: $WORK/flip.evec" \
  "projectiononly: YES" "hiprec: YES"
"$FEAT" -p "$WORK/flip.par" > "$WORK/flip.log" 2>&1
grep -E "recoded|SNP check" "$WORK/flip.log" || true
"$PY" "$TESTS/compare_evec.py" "$WORK/baseline.evec" "$WORK/flip.evec" --tol "$TOL" --npc "$NPC" \
  | grep -E "maximum abs|RESULT" || rc=1
if "$PY" "$TESTS/compare_evec.py" "$WORK/baseline.evec" "$WORK/flip.evec" --tol "$TOL" --npc "$NPC" >/dev/null 2>&1; then :; else rc=1; fi

# ---- validation guard: strict mode must FAIL when a model SNP is absent.
echo
echo "---- validation guard (strict mode rejects incompatible SNP set) ----"
"$PY" - "$DATA/test.snp" "$DATA/test.geno" "$DATA/test.ind" "$WORK/drop" <<'EOF'
import sys
snp,geno,ind,out=sys.argv[1:5]
s=open(snp).read().splitlines(); g=open(geno).read().splitlines()
open(out+".snp","w").write("\n".join(s[50:])+"\n")   # drop 50 model SNPs
open(out+".geno","w").write("\n".join(g[50:])+"\n")
open(out+".ind","w").write(open(ind).read())
EOF
gen_par "$WORK/drop.par" \
  "genotypename: $WORK/drop.geno" \
  "snpname: $WORK/drop.snp" \
  "indivname: $WORK/drop.ind" \
  "modelname: $WORK/reference.smartpca.model" \
  "projectionoutname: $WORK/drop.evec" \
  "projectiononly: YES"
if "$FEAT" -p "$WORK/drop.par" > "$WORK/drop.log" 2>&1; then
  echo "strict-mode guard: FAIL (should have refused to project)"; rc=1
else
  echo "strict-mode guard: PASS (refused, exit nonzero)"
  grep -iE "mismatch|absent" "$WORK/drop.log" | head -1 || true
fi

echo
if [ "$rc" -eq 0 ]; then
  echo "ALL TESTS PASSED"
else
  echo "SOME TESTS FAILED"
fi
exit "$rc"
