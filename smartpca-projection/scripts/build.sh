#!/usr/bin/env bash
#
# build.sh -- build smartpca (EIGENSOFT v18140) with the persistent-projection
# feature.
#
# Strategy, in order of preference:
#   1. If OpenBLAS/LAPACK (and, ideally, GSL) are available -- e.g. inside the
#      conda environment created from environment.yml -- build the full, stock
#      EIGENSOFT smartpca via its own Makefile.  This is the production path and
#      matches upstream numerics exactly.
#   2. If GSL is missing, build without the fast-PCA (fastmode) objects, using a
#      tiny stub for kjg_fpca().  The classic lsqproject PCA + projection path
#      (which this feature targets) is unaffected.
#   3. If no BLAS/LAPACK is available at all, build a local static Reference
#      LAPACK (cloned from github.com/Reference-LAPACK/lapack) and link against
#      it.  This keeps the build self-contained for CI/sandbox use.
#
# Output binary: build/smartpca
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRCROOT="$HERE/upstream"
SRC="$SRCROOT/src"
BUILD="$HERE/build"
mkdir -p "$BUILD"

CC="${CC:-cc}"
FC="${FC:-gfortran}"

echo "== smartpca persistent-projection build =="
echo "   source : $SRC"
echo "   output : $BUILD/smartpca"

# ---------------------------------------------------------------- detect deps
have_gsl=0
have_blas=0
GSL_CFLAGS=""; GSL_LIBS=""
BLAS_LIBS=""

if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists gsl 2>/dev/null; then
  have_gsl=1
  GSL_CFLAGS="$(pkg-config --cflags gsl)"
  GSL_LIBS="$(pkg-config --libs gsl)"
fi

# Prefer OpenBLAS (bundles LAPACK). Fall back to separate lapack+blas.
if [ -n "${CONDA_PREFIX:-}" ] && ls "$CONDA_PREFIX"/lib/libopenblas* >/dev/null 2>&1; then
  have_blas=1; BLAS_LIBS="-L$CONDA_PREFIX/lib -lopenblas"
elif ls /usr/lib*/libopenblas* /usr/lib/*/libopenblas* >/dev/null 2>&1; then
  have_blas=1; BLAS_LIBS="-lopenblas"
elif ls /usr/lib*/liblapack* /usr/lib/*/liblapack* >/dev/null 2>&1; then
  have_blas=1; BLAS_LIBS="-llapack -lblas"
fi

# ------------------------------------------------ build reference LAPACK (opt)
if [ "$have_blas" -eq 0 ]; then
  echo "-- no system BLAS/LAPACK found; building bundled Reference-LAPACK"
  LAPACK_DIR="$BUILD/lapack"
  if [ ! -f "$LAPACK_DIR/liblapack.a" ]; then
    if [ ! -d "$LAPACK_DIR/.git" ] && [ ! -f "$LAPACK_DIR/Makefile" ]; then
      git clone --depth 1 https://github.com/Reference-LAPACK/lapack.git "$LAPACK_DIR"
    fi
    cp "$LAPACK_DIR/make.inc.example" "$LAPACK_DIR/make.inc"
    make -C "$LAPACK_DIR" -j"$(nproc 2>/dev/null || echo 4)" blaslib lapacklib
  fi
  BLAS_LIBS="$LAPACK_DIR/liblapack.a $LAPACK_DIR/librefblas.a -lgfortran"
  have_blas=1
fi

echo "-- GSL:  $([ $have_gsl -eq 1 ] && echo yes || echo 'no (fastmode disabled via stub)')"
echo "-- BLAS: $BLAS_LIBS"

# ---------------------------------------------------------------- compile
CFLAGS="-O2 -g -Wimplicit -I$SRCROOT/include -I$SRC/nicksrc $GSL_CFLAGS"

# nicklib (needs the include dir added explicitly)
make -C "$SRC/nicksrc" CFLAGS="-c -O3 -g -Wimplicit -I./ -I$SRCROOT/include"

OBJS=""
compile () { $CC $CFLAGS -c -o "$SRC/$1.o" "$SRC/$1.c"; OBJS="$OBJS $SRC/$1.o"; }

compile eigensrc/smartpca
compile eigensrc/smartpca_model        # NEW: persistent projection model I/O
compile eigensrc/eigsubs
compile eigensrc/exclude
compile eigensrc/smartsubs
compile eigensrc/eigx
compile mcio
compile qpsubs
compile admutils
compile egsubs
compile regsubs
compile gval

if [ "$have_gsl" -eq 1 ]; then
  # real fast-PCA path
  $CC $CFLAGS -c -o "$SRC/ksrc/kjg_fpca.o" "$SRC/ksrc/kjg_fpca.c"
  $CC $CFLAGS -c -o "$SRC/ksrc/kjg_gsl.o"  "$SRC/ksrc/kjg_gsl.c"
  KOBJS="$SRC/ksrc/kjg_fpca.o $SRC/ksrc/kjg_gsl.o"
else
  # fastmode disabled: link a stub that aborts if fastmode is ever requested
  $CC -O2 -g -I"$SRCROOT/include" -c -o "$SRC/ksrc/kjg_fpca_stub.o" "$SRC/ksrc/kjg_fpca_stub.c"
  KOBJS="$SRC/ksrc/kjg_fpca_stub.o"
  GSL_LIBS=""
fi

# ---------------------------------------------------------------- link
$CC -O2 -g -o "$BUILD/smartpca" \
  $OBJS $KOBJS \
  "$SRC/nicksrc/libnick.a" \
  $BLAS_LIBS $GSL_LIBS \
  -lm -lpthread

echo "== built $BUILD/smartpca =="
"$BUILD/smartpca" -v 2>&1 | sed -n '1,3p' || true
