#!/usr/bin/env python3
"""
Generate a synthetic EIGENSTRAT dataset for smartpca projection equivalence
tests.  Pure-Python (no numpy) so it runs in a minimal environment.

Creates a reference panel of several populations drawn from an admixture-like
model with structured allele frequencies, plus a set of "test" individuals to
be projected -- provided at three missingness levels (complete, moderate, high).

Outputs EIGENSTRAT triples:
  <prefix>.geno  <prefix>.snp  <prefix>.ind

Two datasets are written:
  ref+test combined (for the stock baseline run)   -> <out>/all.{geno,snp,ind}
  reference only                                    -> <out>/ref.{geno,snp,ind}
  test only                                         -> <out>/test.{geno,snp,ind}
and a poplist of the reference populations          -> <out>/ref.pops
"""
import argparse
import os
import random


def make_freqs(nsnps, npops, rng):
    """Ancestral freq per SNP + per-population drifted freqs."""
    anc = [rng.uniform(0.1, 0.9) for _ in range(nsnps)]
    pops = []
    for _ in range(npops):
        fst = rng.uniform(0.01, 0.08)
        f = []
        for a in anc:
            # Balding-Nichols style drift
            alpha = a * (1 - fst) / fst
            beta = (1 - a) * (1 - fst) / fst
            f.append(min(0.99, max(0.01, rng.betavariate(alpha, beta))))
        pops.append(f)
    return pops


def draw_geno(freq, rng):
    """Diploid genotype dosage (count of variant allele) 0/1/2."""
    return (1 if rng.random() < freq else 0) + (1 if rng.random() < freq else 0)


def write_eigenstrat(prefix, geno_rows, snp_meta, ind_meta):
    # geno_rows: list (per SNP) of list of dosages (per individual), -1 = missing
    with open(prefix + ".geno", "w") as g:
        for row in geno_rows:
            g.write("".join(("9" if x < 0 else str(x)) for x in row) + "\n")
    with open(prefix + ".snp", "w") as s:
        for (sid, chrom, gpos, ppos, a0, a1) in snp_meta:
            s.write("%20s %2d %12.6f %12d %s %s\n"
                    % (sid, chrom, gpos, ppos, a0, a1))
    with open(prefix + ".ind", "w") as i:
        for (iid, sex, pop) in ind_meta:
            i.write("%20s %s %s\n" % (iid, sex, pop))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--nsnps", type=int, default=4000)
    ap.add_argument("--npops", type=int, default=5)
    ap.add_argument("--per-pop", type=int, default=20)
    ap.add_argument("--ntest", type=int, default=12)
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.out, exist_ok=True)

    nsnps = args.nsnps
    npops = args.npops
    popfreqs = make_freqs(nsnps, npops, rng)

    # SNP metadata: spread across 22 autosomes
    snp_meta = []
    for k in range(nsnps):
        chrom = (k % 22) + 1
        ppos = 100000 + k * 5000
        gpos = ppos / 1.0e8
        snp_meta.append(("rs%06d" % k, chrom, gpos, ppos, "A", "G"))

    # Reference individuals
    ref_ind = []
    ref_pop_of = []
    for p in range(npops):
        for j in range(args.per_pop):
            ref_ind.append(("REF_P%d_%03d" % (p, j), "U", "Pop%d" % p))
            ref_pop_of.append(p)
    nref = len(ref_ind)

    # Test individuals: assign each to a random reference population ancestry.
    test_pop_of = [rng.randrange(npops) for _ in range(args.ntest)]
    # missingness plan: first third complete, second third moderate (~30%),
    # last third high (~80%)
    def miss_rate(idx):
        third = args.ntest / 3.0
        if idx < third:
            return 0.0
        if idx < 2 * third:
            return 0.30
        return 0.80

    test_ind = []
    for t in range(args.ntest):
        mr = miss_rate(t)
        tag = "COMPLETE" if mr == 0 else ("MOD%02d" % int(mr * 100)) \
            if mr < 0.5 else ("HIGH%02d" % int(mr * 100))
        test_ind.append(("TEST_%02d_%s" % (t, tag), "U", "Test"))

    # Build genotype matrices (per SNP row)
    ref_geno = []
    test_geno = []
    for k in range(nsnps):
        rrow = []
        for i in range(nref):
            rrow.append(draw_geno(popfreqs[ref_pop_of[i]][k], rng))
        ref_geno.append(rrow)

        trow = []
        for t in range(args.ntest):
            g = draw_geno(popfreqs[test_pop_of[t]][k], rng)
            if rng.random() < miss_rate(t):
                g = -1
            trow.append(g)
        test_geno.append(trow)

    # Reference-only dataset
    write_eigenstrat(os.path.join(args.out, "ref"),
                     ref_geno, snp_meta, ref_ind)
    # Test-only dataset
    write_eigenstrat(os.path.join(args.out, "test"),
                     test_geno, snp_meta, test_ind)
    # Combined (reference + test) dataset for the stock baseline run
    all_geno = [ref_geno[k] + test_geno[k] for k in range(nsnps)]
    write_eigenstrat(os.path.join(args.out, "all"),
                     all_geno, snp_meta, ref_ind + test_ind)

    # poplist = the reference populations only (these build the PCA)
    with open(os.path.join(args.out, "ref.pops"), "w") as f:
        for p in range(npops):
            f.write("Pop%d\n" % p)

    print("wrote synthetic data to %s: nsnps=%d nref=%d ntest=%d npops=%d"
          % (args.out, nsnps, nref, args.ntest, npops))


if __name__ == "__main__":
    main()
